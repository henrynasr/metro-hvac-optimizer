"""
sobol_B.py — Sobol sensitivity analysis on the water regime.

The water regime params (T_CW_SUPPLY_MIN/MAX, T_HW_SUPPLY_MIN/MAX,
T_BLOW_COOL, T_BLOW_HEAT, FRAC_RETURN_AIR, DT_WATER_COOL, DT_WATER_HEAT,
T_CW_EXT_LOW, T_CW_EXT_HYST, T_HW_EXT_HIGH, T_HW_EXT_HYST) do NOT feed
dT_dt — they don't affect T_in at all. They only affect Q_water, which is
a post-hoc energy cost metric computed from T_in and n_people.

Strategy:
  1. Run the ODE once at baseline → get T_in_base and n_people_base.
  2. For each Saltelli row: recompute Q_water_cool and Q_water_heat
     post-hoc using the swept params. No ODE solve needed.
     Runtime: ~milliseconds per row instead of ~10ms.

Outputs:
  Y1 — peak Q_water_cool [m³/s]   (chiller sizing)
  Y2 — cumulative Q_water_cool [m³] over the week  (energy cost proxy)
  Y3 — peak Q_water_heat [m³/s]   (boiler sizing)

Ordering constraint:
  T_CW_EXT_LOW < T_CW_EXT_HYST  → swap if violated
  T_HW_EXT_HYST < T_HW_EXT_HIGH → swap if violated
"""

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
from regulation import dT_dt, airflow_total, T_setpoint
from constants import (
    RHO_AIR_KG_M3, CP_AIR_J_KG_K,
    RHO_GLYCOL_KG_M3, CP_GLYCOL_J_KG_K,
    UA_FACADE_W_K, UA_SOIL_W_K, C_TOTAL_J_K,
    T_HW_EXT_LOW_C, T_CW_EXT_HIGH_C,
    AIRFLOW_OVERPRESSURE_M3H, PEOPLE_PEAK,
)

RHO_CP = RHO_AIR_KG_M3 * CP_AIR_J_KG_K


# =============================================================================
# SECTION 1 — Baseline ODE solve (once)
# =============================================================================

def run_baseline_ode(t_array, T_ext_array, Q_int_array, n_people_array,
                     dates, T0=23.0):
    """
    Solve the thermal ODE once at baseline constants.
    Returns T_in_array over the simulation window.
    """
    water_state = {"heating": False, "cooling": False}
    sol = solve_ivp(
        dT_dt,
        t_span  = (t_array[0], t_array[-1]),
        y0      = [T0],
        t_eval  = t_array,
        args    = (t_array, T_ext_array, Q_int_array, n_people_array,
                   dates, water_state),
        method  = "RK45",
        max_step= 3600.0,
    )
    return sol.y[0]


# =============================================================================
# SECTION 2 — Parametric Q_water computation (no ODE)
# =============================================================================

def eval_Q_water(row, param_names,
                 T_in_array, T_ext_array, n_people_array, dates):
    """
    Recompute Q_water_cool and Q_water_heat for one parameter set,
    given a fixed T_in trajectory (baseline ODE output).

    No ODE solve — post-hoc only.

    Returns
    -------
    (peak_Qwc, cumul_Qwc, peak_Qwh) : (float, float, float)
    """
    p = dict(zip(param_names, row))

    # Ordering constraint guard
    if p["T_HW_EXT_HYST"] >= p["T_HW_EXT_HIGH"]:
        p["T_HW_EXT_HYST"], p["T_HW_EXT_HIGH"] = p["T_HW_EXT_HIGH"], p["T_HW_EXT_HYST"]
    if p["T_CW_EXT_LOW"] >= p["T_CW_EXT_HYST"]:
        p["T_CW_EXT_LOW"], p["T_CW_EXT_HYST"] = p["T_CW_EXT_HYST"], p["T_CW_EXT_LOW"]

    n            = len(T_in_array)
    Q_water_cool = np.zeros(n)
    Q_water_heat = np.zeros(n)

    heating_active = False
    cooling_active = False

    for i, (T_in, T_ext, n_ppl) in enumerate(
            zip(T_in_array, T_ext_array, n_people_array)):

        T_set     = T_setpoint(T_ext)          # baseline setpoint law (not swept here)
        Q_air_m3s = airflow_total(n_ppl) / 3600.0
        T_mix     = p["FRAC_RETURN_AIR"] * T_in + (1.0 - p["FRAC_RETURN_AIR"]) * T_ext

        # --- Hot water hysteresis ---
        if heating_active:
            if T_ext >= p["T_HW_EXT_HIGH"]:
                heating_active = False
        else:
            if T_ext <= p["T_HW_EXT_HYST"]:
                heating_active = True

        # --- Cold water hysteresis ---
        if cooling_active:
            if T_ext <= p["T_CW_EXT_LOW"]:
                cooling_active = False
        else:
            if T_ext >= p["T_CW_EXT_HYST"]:
                cooling_active = True

        # --- Hot water supply temp (linear interp, swept min/max) ---
        if heating_active:
            T_hw = float(np.interp(
                T_ext,
                [T_HW_EXT_LOW_C,        p["T_HW_EXT_HIGH"]],
                [p["T_HW_SUPPLY_MAX"],   p["T_HW_SUPPLY_MIN"]],
            ))
        else:
            T_hw = None

        # --- Cold water supply temp (linear interp, swept min/max) ---
        if cooling_active:
            T_cw = float(np.interp(
                T_ext,
                [p["T_CW_EXT_LOW"],      T_CW_EXT_HIGH_C],
                [p["T_CW_SUPPLY_MAX"],   p["T_CW_SUPPLY_MIN"]],
            ))
        else:
            T_cw = None

        # --- Q_water ---
        if not np.isnan(T_set) and T_set > T_in and heating_active:
            dT_air = p["T_BLOW_HEAT"] - T_mix
            if dT_air > 0:
                Q_water_heat[i] = (
                    Q_air_m3s * RHO_CP * dT_air
                    / (RHO_GLYCOL_KG_M3 * CP_GLYCOL_J_KG_K * p["DT_WATER_HEAT"])
                )

        if not np.isnan(T_set) and T_set <= T_in and cooling_active:
            dT_air = T_mix - p["T_BLOW_COOL"]
            if dT_air > 0:
                Q_water_cool[i] = (
                    Q_air_m3s * RHO_CP * dT_air
                    / (RHO_GLYCOL_KG_M3 * CP_GLYCOL_J_KG_K * p["DT_WATER_COOL"])
                )

    # Outputs
    peak_Qwc  = float(Q_water_cool.max())
    cumul_Qwc = float(Q_water_cool.sum() * 3600.0)   # m³ (each step = 1h = 3600s)
    peak_Qwh  = float(Q_water_heat.max())

    return peak_Qwc, cumul_Qwc, peak_Qwh


# =============================================================================
# SECTION 3 — Problem definition
# =============================================================================

PROBLEM_B = {
    "num_vars": 12,
    "names": [
        "T_HW_SUPPLY_MAX",   # hot water supply at T_ext = -7°C
        "T_HW_SUPPLY_MIN",   # hot water supply at T_ext = T_HW_EXT_HIGH
        "T_HW_EXT_HIGH",     # hot water shutoff T_ext
        "T_HW_EXT_HYST",     # hot water restart T_ext  (must be < T_HW_EXT_HIGH)
        "T_CW_SUPPLY_MAX",   # cold water supply at T_ext = T_CW_EXT_LOW
        "T_CW_SUPPLY_MIN",   # cold water supply at T_ext = 31°C
        "T_CW_EXT_LOW",      # cold water onset T_ext   (must be < T_CW_EXT_HYST)
        "T_CW_EXT_HYST",     # cold water restart T_ext
        "T_BLOW_COOL",       # AHU supply temp, cooling mode
        "T_BLOW_HEAT",       # AHU supply temp, heating mode
        "FRAC_RETURN_AIR",   # return air fraction (affects T_mix → Q_water)
        "DT_WATER_COOL",     # supply−return ΔT, cold circuit
        "DT_WATER_HEAT",     # supply−return ΔT, hot circuit
    ],
    "bounds": [
        [45.0,  55.0],    # T_HW_SUPPLY_MAX   [°C]
        [28.0,  40.0],    # T_HW_SUPPLY_MIN   [°C]
        [10.0,  15.0],    # T_HW_EXT_HIGH     [°C]
        [8.0,   12.0],    # T_HW_EXT_HYST     [°C]
        [10.0,  16.0],    # T_CW_SUPPLY_MAX   [°C]
        [6.0,   12.0],    # T_CW_SUPPLY_MIN   [°C]
        [24.0,  28.0],    # T_CW_EXT_LOW      [°C]
        [26.0,  30.0],    # T_CW_EXT_HYST     [°C]
        [13.0,  17.0],    # T_BLOW_COOL       [°C]
        [28.0,  35.0],    # T_BLOW_HEAT       [°C]
        [0.60,  0.80],    # FRAC_RETURN_AIR   [-]
        [4.0,   8.0],     # DT_WATER_COOL     [K]
        [3.0,   7.0],     # DT_WATER_HEAT     [K]
    ],
}
PROBLEM_B["num_vars"] = len(PROBLEM_B["names"])


# =============================================================================
# SECTION 4 — Plot
# =============================================================================

def plot_sobol(Si, names, metric_name, filename):
    x     = np.arange(len(names))
    width = 0.35
    fig, ax = plt.subplots(figsize=(max(10, len(names) * 0.8), 6))
    ax.bar(x - width/2, Si["S1"], width, yerr=Si["S1_conf"],
           label="S1 (first-order)", capsize=4, color="steelblue")
    ax.bar(x + width/2, Si["ST"], width, yerr=Si["ST_conf"],
           label="ST (total effect)", capsize=4, color="firebrick")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=9)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.legend()
    style_axes(ax, title=f"Sobol B — {metric_name}", ylabel="Sobol index")
    plt.tight_layout()
    plt.savefig(filename, dpi=120)
    plt.close()
    print(f"  Saved {filename}")


# =============================================================================
# SECTION 5 — Main
# =============================================================================

if __name__ == "__main__":

    print("Loading data...")
    df      = load_data("data/raw/paris_weather.csv")
    df_week = df.loc["2024-08-01":"2024-08-07"]

    t_array     = np.arange(len(df_week)) * 3600.0
    T_ext_array = df_week["temperature_2m"].values
    dates       = df_week.index

    profiles           = load_profiles("data/raw/Defense_Occupation_Normalised.xlsx")
    Q_int_base, n_ppl_base = build_Q_array(dates, profiles)

    # --- Step 1: baseline ODE solve ---
    print("Running baseline ODE...")
    T_in_base = run_baseline_ode(t_array, T_ext_array, Q_int_base,
                                  n_ppl_base, dates)
    print(f"  Baseline T_in: min={T_in_base.min():.1f}°C  "
          f"max={T_in_base.max():.1f}°C")

    # --- Step 2: Saltelli sampling ---
    N = 512
    samples_B = sobol_sample.sample(PROBLEM_B, N)
    print(f"\nSobol B — {PROBLEM_B['num_vars']} params")
    print(f"  Sample shape: {samples_B.shape}  ({samples_B.shape[0]} evaluations)")

    N_B      = samples_B.shape[0]
    Y_peak_cool  = np.full(N_B, np.nan)
    Y_cumul_cool = np.full(N_B, np.nan)
    Y_peak_heat  = np.full(N_B, np.nan)

    # --- Step 3: main loop (no ODE — fast) ---
    t0 = time.perf_counter()
    for i, row in enumerate(tqdm(samples_B, desc="Sobol B")):
        Y_peak_cool[i], Y_cumul_cool[i], Y_peak_heat[i] = eval_Q_water(
            row, PROBLEM_B["names"],
            T_in_base, T_ext_array, n_ppl_base, dates,
        )
    elapsed = time.perf_counter() - t0
    print(f"\n  Done in {elapsed:.1f} s  ({elapsed/N_B*1000:.2f} ms/eval)")
    print(f"  Peak Q_water_cool: {np.nanmin(Y_peak_cool)*1e6:.1f} – "
          f"{np.nanmax(Y_peak_cool)*1e6:.1f} cm³/s")
    print(f"  Cumul Q_water_cool: {np.nanmin(Y_cumul_cool):.3f} – "
          f"{np.nanmax(Y_cumul_cool):.3f} m³")

    np.save("data/processed/sobol_B_peak_cool.npy",  Y_peak_cool)
    np.save("data/processed/sobol_B_cumul_cool.npy", Y_cumul_cool)
    np.save("data/processed/sobol_B_peak_heat.npy",  Y_peak_heat)

    # --- Step 4: analysis and plots ---
    Si_peak_cool  = sobol_analyze.analyze(PROBLEM_B, Y_peak_cool,  print_to_console=False)
    Si_cumul_cool = sobol_analyze.analyze(PROBLEM_B, Y_cumul_cool, print_to_console=False)
    Si_peak_heat  = sobol_analyze.analyze(PROBLEM_B, Y_peak_heat,  print_to_console=False)

    plot_sobol(Si_peak_cool,  PROBLEM_B["names"],
               "Peak Q_water_cool [m³/s]",    "images/sobol_B_peak_cool.png")
    plot_sobol(Si_cumul_cool, PROBLEM_B["names"],
               "Cumulative Q_water_cool [m³]", "images/sobol_B_cumul_cool.png")
    plot_sobol(Si_peak_heat,  PROBLEM_B["names"],
               "Peak Q_water_heat [m³/s]",    "images/sobol_B_peak_heat.png")

    # --- Step 5: print top drivers ---
    print(f"\n{'='*50}")
    print("TOP DRIVERS — Peak Q_water_cool (ST index)")
    ranked = sorted(zip(PROBLEM_B["names"], Si_peak_cool["ST"]),
                    key=lambda x: x[1], reverse=True)
    for name, st in ranked:
        print(f"  {name:<22} ST = {st:.3f}")
