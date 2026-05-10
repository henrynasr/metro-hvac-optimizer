# =============================================================================
# sobol_A.py — Sobol GSA on the thermal model
# Run A: 26 params, full screen, N=128
# Run C: 5 survivors from A, geometry fixed, N=1024
# =============================================================================

import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from tqdm import tqdm
from SALib.sample import sobol as sobol_sample
from SALib.analyze import sobol as sobol_analyze

from utils import load_data, style_axes
from occupancy import load_profiles, build_Q_array, _day_type
from constants import (
    RHO_AIR_KG_M3, CP_AIR_J_KG_K, RHO_CONC_KG_M3, CP_CONC_J_KG_K,
    RHO_GLYCOL_KG_M3, CP_GLYCOL_J_KG_K,
    PLATFORM_LENGTH_M, AIRFLOW_PER_PERSON_M3H,
    T_HW_EXT_LOW_C,
    U_FACADE_W_M2K, U_SOIL_W_M2K, D_CONC_EFF_M,
    T_TUN_OFFSET_DAY_C, T_TUN_OFFSET_NIGHT_C, T_SOIL_C, ETA_INF,
    HEADWAY_PEAK_S, HEADWAY_OFFPEAK_S,
    PEOPLE_PEAK, WATTS_SENSIBLE_PP, BASELINE_W,
    AIRFLOW_OVERPRESSURE_M3H, T_BLOW_COOL_C, T_BLOW_HEAT_C, FRAC_RETURN_AIR,
    T_ANTIFREEZE_C, T_DEAD_LOW_C, T_DEAD_HIGH_C, T_COOL_DELTA_C,
    T_HW_EXT_HIGH_C, T_HW_EXT_HYST_C, T_CW_EXT_LOW_C, T_CW_EXT_HYST_C,
    UA_TUN_WALL_W_K, Q_STAIR_M3S,
)

RHO_CP = RHO_AIR_KG_M3 * CP_AIR_J_KG_K   # 1206 J/(m³·K)


# -----------------------------------------------------------------------------
# 1. DEFAULTS — baseline fallback for params not swept in a given run
# -----------------------------------------------------------------------------

DEFAULTS = {
    "W":                    4.0,
    "H":                    4.2,
    "H_F":                  2.8,
    "U_FACADE":             U_FACADE_W_M2K,
    "U_SOIL":               U_SOIL_W_M2K,
    "D_CONC_EFF":           D_CONC_EFF_M,
    "T_TUN_OFFSET_DAY":     T_TUN_OFFSET_DAY_C,
    "T_TUN_OFFSET_NIGHT":   T_TUN_OFFSET_NIGHT_C,
    "T_SOIL":               T_SOIL_C,
    "ETA_INF":              ETA_INF,
    "HW_PEAK":              HEADWAY_PEAK_S,
    "HW_OFFPEAK":           HEADWAY_OFFPEAK_S,
    "PEOPLE_PEAK":          float(PEOPLE_PEAK),
    "WATTS_PP":             WATTS_SENSIBLE_PP,
    "BASELINE_W":           BASELINE_W,
    "AIRFLOW_OVERPRESSURE": AIRFLOW_OVERPRESSURE_M3H,
    "T_BLOW_COOL":          T_BLOW_COOL_C,
    "T_BLOW_HEAT":          T_BLOW_HEAT_C,
    "FRAC_RETURN_AIR":      FRAC_RETURN_AIR,
    "T_ANTIFREEZE":         T_ANTIFREEZE_C,
    "T_DEAD_LOW":           T_DEAD_LOW_C,
    "T_DEAD_HIGH":          T_DEAD_HIGH_C,
    "T_COOL_DELTA":         T_COOL_DELTA_C,
    "T_HW_EXT_HIGH":        T_HW_EXT_HIGH_C,
    "T_HW_EXT_HYST":        T_HW_EXT_HYST_C,
    "T_CW_EXT_LOW":         T_CW_EXT_LOW_C,
    "T_CW_EXT_HYST":        T_CW_EXT_HYST_C,
}


# -----------------------------------------------------------------------------
# 2. PARAMETRIC HELPERS
# -----------------------------------------------------------------------------

def _v_inf(hour, day_type, eta, V_platform, hw_peak, hw_offpeak):
    WEEKDAY = {**{h:"night" for h in range(5)}, 5:"offpeak",
               **{h:"peak" for h in range(6,9)}, **{h:"offpeak" for h in range(9,17)},
               **{h:"peak" for h in range(17,20)}, **{h:"offpeak" for h in range(20,23)}, 23:"night"}
    WKD     = {h: ("offpeak" if 6 <= h <= 22 else "night") for h in range(24)}
    regime  = WKD[hour] if day_type == "WKD" else WEEKDAY[hour]
    V_cycle = eta * V_platform
    if regime == "night":  return 0.0
    if regime == "peak":   return V_cycle / hw_peak
    return V_cycle / hw_offpeak


def _T_set(T_ext, T_antifreeze, T_dead_low, T_dead_high, T_cool_delta):
    """Two-zone setpoint law: heat to 21°C / cool to T_ext-delta."""
    if T_ext < T_dead_low:   return max(21.0, T_antifreeze)
    if T_ext <= T_dead_high: return np.nan
    return T_ext - T_cool_delta


def _airflow(n_people, airflow_overpressure, people_peak):
    Q_min = airflow_overpressure
    Q_max = (airflow_overpressure + people_peak * AIRFLOW_PER_PERSON_M3H) * 1.10
    Q     = airflow_overpressure + min(n_people, people_peak) * AIRFLOW_PER_PERSON_M3H
    return float(np.clip(Q, Q_min, Q_max))


def _ode(t, T, t_array, T_ext_array, Q_int_array, n_people_array, dates,
         UA_facade, UA_tun_wall, UA_soil, C_total,
         T_tun_offset_day, T_tun_offset_night, T_soil,
         eta, V_platform, hw_peak, hw_offpeak,
         people_peak, airflow_overpressure,
         T_blow_heat, T_blow_cool,
         T_antifreeze, T_dead_low, T_dead_high, T_cool_delta,
         water_state,
         T_hw_ext_high, T_hw_ext_hyst, T_cw_ext_low, T_cw_ext_hyst,
         Q_stair_m3s):
    """Parametric ODE — all constants explicit, no module imports."""
    T_in  = float(T[0])
    T_ext = float(np.interp(t, t_array, T_ext_array))
    Q_int = float(np.interp(t, t_array, Q_int_array))
    n_ppl = float(np.interp(t, t_array, n_people_array))

    idx   = int(np.clip(np.searchsorted(t_array, t), 0, len(dates) - 1))
    ts    = dates[idx];  dtype = _day_type(ts);  hour = ts.hour

    offset     = T_tun_offset_day if 5 <= hour <= 23 else T_tun_offset_night
    T_tun      = T_ext + offset
    V_inf      = _v_inf(hour, dtype, eta, V_platform, hw_peak, hw_offpeak)
    UA_eff_tun = UA_facade + UA_tun_wall + RHO_CP * V_inf

    Q_facade = UA_eff_tun * (T_tun  - T_in)
    Q_soil   = UA_soil    * (T_soil - T_in)
    Q_stair  = RHO_CP * Q_stair_m3s * (T_ext - T_in)

    Tset      = _T_set(T_ext, T_antifreeze, T_dead_low, T_dead_high, T_cool_delta)
    Q_air_m3s = _airflow(n_ppl, airflow_overpressure, people_peak) / 3600.0

    hon = water_state["heating"];  con = water_state["cooling"]
    if hon:
        if T_ext >= T_hw_ext_high: hon = False
    else:
        if T_ext <= T_hw_ext_hyst: hon = True
    if con:
        if T_ext <= T_cw_ext_low:  con = False
    else:
        if T_ext >= T_cw_ext_hyst: con = True
    water_state["heating"] = hon;  water_state["cooling"] = con

    if np.isnan(Tset):
        Q_hvac = 0.0
    elif Tset > T_in:
        Q_hvac = RHO_CP * Q_air_m3s * (T_in - T_blow_heat)
        if not hon: Q_hvac = 0.0
    else:
        Q_hvac = RHO_CP * Q_air_m3s * (T_in - T_blow_cool)
        if not con: Q_hvac = 0.0

    return [(Q_facade + Q_soil + Q_int + Q_stair - Q_hvac) / C_total]


# -----------------------------------------------------------------------------
# 3. GEOMETRY HELPERS
# -----------------------------------------------------------------------------

def _geometry(W, H, H_f):
    L = PLATFORM_LENGTH_M
    return {
        "A_facade":   L * H_f,
        "A_soil":     L * H + W * H + 2 * L * W,   # one end wall (stair side excluded)
        "V_platform": L * W * H,
    }

def _capacitance(A_soil, d_eff, V_platform):
    return (RHO_AIR_KG_M3  * CP_AIR_J_KG_K  * V_platform
          + RHO_CONC_KG_M3 * CP_CONC_J_KG_K * A_soil * d_eff)

_BASELINE_GEO = _geometry(4.0, 4.2, 2.8)


# -----------------------------------------------------------------------------
# 4. EVAL MODEL
# -----------------------------------------------------------------------------

def eval_model(row, param_names, t_array, T_ext_array,
               dates, profiles, fixed_geo=None, T0=23.0):
    """One ODE solve per Saltelli row. fixed_geo=None → sweep geometry (run A)."""
    p = {**DEFAULTS, **dict(zip(param_names, row))}

    # Ordering constraint guards
    if p["T_HW_EXT_HYST"] >= p["T_HW_EXT_HIGH"]:
        p["T_HW_EXT_HYST"], p["T_HW_EXT_HIGH"] = p["T_HW_EXT_HIGH"], p["T_HW_EXT_HYST"]
    if p["T_CW_EXT_LOW"] >= p["T_CW_EXT_HYST"]:
        p["T_CW_EXT_LOW"], p["T_CW_EXT_HYST"] = p["T_CW_EXT_HYST"], p["T_CW_EXT_LOW"]

    geo     = fixed_geo if fixed_geo is not None else _geometry(p["W"], p["H"], p["H_F"])
    UA_f    = p["U_FACADE"] * geo["A_facade"]
    UA_s    = p["U_SOIL"]   * geo["A_soil"]
    C_total = _capacitance(geo["A_soil"], p["D_CONC_EFF"], geo["V_platform"])

    Q_int_sw, n_ppl_sw = build_Q_array(
        dates, profiles,
        watts_per_person = p["WATTS_PP"],
        baseline_w       = p["BASELINE_W"],
        people_peak      = int(p["PEOPLE_PEAK"]),
    )

    water_state = {"heating": False, "cooling": False}

    sol = solve_ivp(
        _ode,
        t_span   = (t_array[0], t_array[-1]),
        y0       = [T0],
        t_eval   = t_array,
        method   = "RK45",
        max_step = 3600.0,
        args     = (
            t_array, T_ext_array, Q_int_sw, n_ppl_sw, dates,
            UA_f, UA_TUN_WALL_W_K, UA_s, C_total,
            p["T_TUN_OFFSET_DAY"], p["T_TUN_OFFSET_NIGHT"], p["T_SOIL"],
            p["ETA_INF"], geo["V_platform"], p["HW_PEAK"], p["HW_OFFPEAK"],
            p["PEOPLE_PEAK"], p["AIRFLOW_OVERPRESSURE"],
            p["T_BLOW_HEAT"], p["T_BLOW_COOL"],
            p["T_ANTIFREEZE"], p["T_DEAD_LOW"], p["T_DEAD_HIGH"], p["T_COOL_DELTA"],
            water_state,
            p["T_HW_EXT_HIGH"], p["T_HW_EXT_HYST"],
            p["T_CW_EXT_LOW"],  p["T_CW_EXT_HYST"],
            Q_STAIR_M3S,
        ),
    )

    T_in = sol.y[0]
    return float(T_in.max()), float((T_in > 26).sum() / len(T_in) * 100)


# -----------------------------------------------------------------------------
# 5. PROBLEM DEFINITIONS
# -----------------------------------------------------------------------------

PROBLEM_A = {
    "names": [
        "W", "H", "H_F",
        "U_FACADE", "U_SOIL", "D_CONC_EFF",
        "T_TUN_OFFSET_DAY", "T_TUN_OFFSET_NIGHT", "T_SOIL",
        "ETA_INF", "HW_PEAK", "HW_OFFPEAK",
        "PEOPLE_PEAK", "WATTS_PP", "BASELINE_W",
        "AIRFLOW_OVERPRESSURE", "T_BLOW_COOL", "T_BLOW_HEAT", "FRAC_RETURN_AIR",
        "T_ANTIFREEZE", "T_DEAD_LOW", "T_DEAD_HIGH", "T_COOL_DELTA",
        "T_HW_EXT_HIGH", "T_HW_EXT_HYST", "T_CW_EXT_LOW", "T_CW_EXT_HYST",
    ],
    "bounds": [
        [3.5,   5.0],    # W
        [3.8,   5.0],    # H
        [2.5,   3.2],    # H_F
        [5.5,   8.5],    # U_FACADE
        [0.35,  0.59],   # U_SOIL
        [0.05,  0.20],   # D_CONC_EFF
        [5.0,   15.0],   # T_TUN_OFFSET_DAY
        [3.0,   8.0],    # T_TUN_OFFSET_NIGHT
        [13.0,  17.0],   # T_SOIL
        [0.10,  0.25],   # ETA_INF
        [90.0,  120.0],  # HW_PEAK
        [180.0, 300.0],  # HW_OFFPEAK
        [150.0, 350.0],  # PEOPLE_PEAK
        [70.0,  90.0],   # WATTS_PP
        [300.,  8000.],  # BASELINE_W
        [1500., 4000.],  # AIRFLOW_OVERPRESSURE
        [13.0,  17.0],   # T_BLOW_COOL
        [28.0,  35.0],   # T_BLOW_HEAT
        [0.60,  0.80],   # FRAC_RETURN_AIR
        [3.0,   7.0],    # T_ANTIFREEZE
        [12.0,  18.0],   # T_DEAD_LOW
        [18.0,  22.0],   # T_DEAD_HIGH
        [4.0,   8.0],    # T_COOL_DELTA
        [10.0,  15.0],   # T_HW_EXT_HIGH
        [8.0,   12.0],   # T_HW_EXT_HYST
        [24.0,  28.0],   # T_CW_EXT_LOW
        [26.0,  30.0],   # T_CW_EXT_HYST
    ],
}
PROBLEM_A["num_vars"] = len(PROBLEM_A["names"])

PROBLEM_C = {
    "names": ["D_CONC_EFF", "T_TUN_OFFSET_DAY", "PEOPLE_PEAK", "BASELINE_W", "AIRFLOW_OVERPRESSURE"],
    "bounds": [
        [0.05,  0.20],
        [5.0,   15.0],
        [150.0, 350.0],
        [300.,  8000.],
        [1500., 4000.],
    ],
}
PROBLEM_C["num_vars"] = len(PROBLEM_C["names"])


# -----------------------------------------------------------------------------
# 6. PLOT HELPER
# -----------------------------------------------------------------------------

def plot_sobol(Si, names, title, filename):
    x, w = np.arange(len(names)), 0.35
    fig, ax = plt.subplots(figsize=(max(12, len(names) * 0.7), 6))
    ax.bar(x - w/2, Si["S1"], w, yerr=Si["S1_conf"], label="S1", capsize=4, color="steelblue")
    ax.bar(x + w/2, Si["ST"], w, yerr=Si["ST_conf"], label="ST", capsize=4, color="firebrick")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=9)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.legend()
    style_axes(ax, title=title, ylabel="Sobol index")
    plt.tight_layout()
    plt.savefig(filename, dpi=120)
    plt.close()
    print(f"  Saved {filename}")


# -----------------------------------------------------------------------------
# 7. MAIN
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    print("Loading data...")
    df      = load_data("data/raw/paris_weather.csv")
    df_week = df.loc["2024-07-01":"2024-07-07"]

    t_array     = np.arange(len(df_week)) * 3600.0
    T_ext_array = df_week["temperature_2m"].values
    dates       = df_week.index
    profiles    = load_profiles("data/raw/Defense_Occupation_Normalised.xlsx")

    # --- Sobol A ---
    N_A       = 16
    samples_A = sobol_sample.sample(PROBLEM_A, N_A)
    print(f"\nSOBOL A — {PROBLEM_A['num_vars']} params | N={N_A} | {samples_A.shape[0]} ODE solves")

    Y_peak_A = np.full(samples_A.shape[0], np.nan)
    Y_pct_A  = np.full(samples_A.shape[0], np.nan)
    t0 = time.perf_counter()
    for i, row in enumerate(tqdm(samples_A, desc="Sobol A")):
        Y_peak_A[i], Y_pct_A[i] = eval_model(row, PROBLEM_A["names"],
                                               t_array, T_ext_array, dates, profiles)
    print(f"  Done in {time.perf_counter()-t0:.1f}s")

    np.save("data/processed/sobol_A_peak.npy", Y_peak_A)
    np.save("data/processed/sobol_A_pct.npy",  Y_pct_A)

    Si_A_peak = sobol_analyze.analyze(PROBLEM_A, Y_peak_A, print_to_console=False)
    Si_A_pct  = sobol_analyze.analyze(PROBLEM_A, Y_pct_A,  print_to_console=False)
    plot_sobol(Si_A_peak, PROBLEM_A["names"], "Sobol A — Peak T_in",      "images/sobol_A_peak.png")
    plot_sobol(Si_A_pct,  PROBLEM_A["names"], "Sobol A — % hours > 26°C", "images/sobol_A_pct.png")

    # --- Sobol C ---
    N_C       = 512
    samples_C = sobol_sample.sample(PROBLEM_C, N_C)
    print(f"\nSOBOL C — {PROBLEM_C['num_vars']} params | N={N_C} | {samples_C.shape[0]} ODE solves")

    Y_peak_C = np.full(samples_C.shape[0], np.nan)
    Y_pct_C  = np.full(samples_C.shape[0], np.nan)
    t0 = time.perf_counter()
    for i, row in enumerate(tqdm(samples_C, desc="Sobol C")):
        Y_peak_C[i], Y_pct_C[i] = eval_model(row, PROBLEM_C["names"],
                                               t_array, T_ext_array, dates, profiles,
                                               fixed_geo=_BASELINE_GEO)
    print(f"  Done in {time.perf_counter()-t0:.1f}s")

    np.save("data/processed/sobol_C_peak.npy", Y_peak_C)
    np.save("data/processed/sobol_C_pct.npy",  Y_pct_C)

    Si_C_peak = sobol_analyze.analyze(PROBLEM_C, Y_peak_C, print_to_console=False)
    Si_C_pct  = sobol_analyze.analyze(PROBLEM_C, Y_pct_C,  print_to_console=False)
    plot_sobol(Si_C_peak, PROBLEM_C["names"], "Sobol C — Peak T_in (top 5)",      "images/sobol_C_peak.png")
    plot_sobol(Si_C_pct,  PROBLEM_C["names"], "Sobol C — % hours > 26°C (top 5)", "images/sobol_C_pct.png")

    print(f"\nSOBOL C — Peak T_in")
    print(f"  {'Param':<25} {'S1':>7} ± conf     {'ST':>7} ± conf")
    for k, name in enumerate(PROBLEM_C["names"]):
        print(f"  {name:<25} {Si_C_peak['S1'][k]:>7.3f} ± {Si_C_peak['S1_conf'][k]:<6.3f}"
              f"   {Si_C_peak['ST'][k]:>7.3f} ± {Si_C_peak['ST_conf'][k]:<6.3f}")