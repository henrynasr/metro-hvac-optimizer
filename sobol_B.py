# =============================================================================
# sobol_B.py — Sobol GSA on water regime parameters (post-hoc, no ODE per row)
# Strategy: run ODE once at baseline → T_in_base, then recompute Q_water
# post-hoc for each Saltelli row. Runtime ~ms/row instead of ~10ms.
# Week: August 2024 (cold circuit fires; July too cold for chiller).
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
from regulation import dT_dt, airflow_total, T_setpoint
from constants import (
    RHO_AIR_KG_M3, CP_AIR_J_KG_K,
    RHO_GLYCOL_KG_M3, CP_GLYCOL_J_KG_K,
    T_HW_EXT_LOW_C, T_CW_EXT_HIGH_C,
)

RHO_CP = RHO_AIR_KG_M3 * CP_AIR_J_KG_K


# -----------------------------------------------------------------------------
# 1. BASELINE ODE SOLVE
# -----------------------------------------------------------------------------

def run_baseline_ode(t_array, T_ext_array, Q_int_array, n_people_array, dates, T0=23.0):
    water_state = {"heating": False, "cooling": False}
    sol = solve_ivp(
        dT_dt,
        t_span=(t_array[0], t_array[-1]), y0=[T0], t_eval=t_array,
        args=(t_array, T_ext_array, Q_int_array, n_people_array, dates, water_state),
        method="RK45", max_step=3600.0,
    )
    return sol.y[0]


# -----------------------------------------------------------------------------
# 2. PARAMETRIC Q_WATER (no ODE)
# -----------------------------------------------------------------------------

def eval_Q_water(row, param_names, T_in_array, T_ext_array, n_people_array, dates):
    """Recompute Q_water for one param set given fixed T_in trajectory."""
    p = dict(zip(param_names, row))

    # Ordering guards
    if p["T_HW_EXT_HYST"] >= p["T_HW_EXT_HIGH"]:
        p["T_HW_EXT_HYST"], p["T_HW_EXT_HIGH"] = p["T_HW_EXT_HIGH"], p["T_HW_EXT_HYST"]
    if p["T_CW_EXT_LOW"] >= p["T_CW_EXT_HYST"]:
        p["T_CW_EXT_LOW"], p["T_CW_EXT_HYST"] = p["T_CW_EXT_HYST"], p["T_CW_EXT_LOW"]

    n = len(T_in_array)
    Q_water_cool = np.zeros(n)
    Q_water_heat = np.zeros(n)
    heating_active = cooling_active = False

    for i, (T_in, T_ext, n_ppl) in enumerate(zip(T_in_array, T_ext_array, n_people_array)):
        T_set     = T_setpoint(T_ext)
        Q_air_m3s = airflow_total(n_ppl) / 3600.0
        T_mix     = p["FRAC_RETURN_AIR"] * T_in + (1.0 - p["FRAC_RETURN_AIR"]) * T_ext

        # Hysteresis
        if heating_active:
            if T_ext >= p["T_HW_EXT_HIGH"]: heating_active = False
        else:
            if T_ext <= p["T_HW_EXT_HYST"]: heating_active = True
        if cooling_active:
            if T_ext <= p["T_CW_EXT_LOW"]:  cooling_active = False
        else:
            if T_ext >= p["T_CW_EXT_HYST"]: cooling_active = True

        if not np.isnan(T_set) and T_set > T_in and heating_active:
            dT_air = p["T_BLOW_HEAT"] - T_mix
            if dT_air > 0:
                Q_water_heat[i] = (Q_air_m3s * RHO_CP * dT_air
                                   / (RHO_GLYCOL_KG_M3 * CP_GLYCOL_J_KG_K * p["DT_WATER_HEAT"]))

        if not np.isnan(T_set) and T_set <= T_in and cooling_active:
            dT_air = T_mix - p["T_BLOW_COOL"]
            if dT_air > 0:
                Q_water_cool[i] = (Q_air_m3s * RHO_CP * dT_air
                                   / (RHO_GLYCOL_KG_M3 * CP_GLYCOL_J_KG_K * p["DT_WATER_COOL"]))

    return float(Q_water_cool.max()), float(Q_water_cool.sum() * 3600.0), float(Q_water_heat.max())


# -----------------------------------------------------------------------------
# 3. PROBLEM DEFINITION
# -----------------------------------------------------------------------------

PROBLEM_B = {
    "names": [
        "T_HW_SUPPLY_MAX", "T_HW_SUPPLY_MIN",
        "T_HW_EXT_HIGH",   "T_HW_EXT_HYST",
        "T_CW_SUPPLY_MAX", "T_CW_SUPPLY_MIN",
        "T_CW_EXT_LOW",    "T_CW_EXT_HYST",
        "T_BLOW_COOL",     "T_BLOW_HEAT",
        "FRAC_RETURN_AIR", "DT_WATER_COOL",   "DT_WATER_HEAT",
    ],
    "bounds": [
        [45.0,  55.0],  # T_HW_SUPPLY_MAX
        [28.0,  40.0],  # T_HW_SUPPLY_MIN
        [10.0,  15.0],  # T_HW_EXT_HIGH
        [8.0,   12.0],  # T_HW_EXT_HYST
        [10.0,  16.0],  # T_CW_SUPPLY_MAX
        [6.0,   12.0],  # T_CW_SUPPLY_MIN
        [24.0,  28.0],  # T_CW_EXT_LOW
        [26.0,  30.0],  # T_CW_EXT_HYST
        [13.0,  17.0],  # T_BLOW_COOL
        [28.0,  35.0],  # T_BLOW_HEAT
        [0.60,  0.80],  # FRAC_RETURN_AIR
        [4.0,   8.0],   # DT_WATER_COOL
        [3.0,   7.0],   # DT_WATER_HEAT
    ],
}
PROBLEM_B["num_vars"] = len(PROBLEM_B["names"])   # 13


# -----------------------------------------------------------------------------
# 4. PLOT HELPER
# -----------------------------------------------------------------------------

def plot_sobol(Si, names, metric_name, filename):
    x, w = np.arange(len(names)), 0.35
    fig, ax = plt.subplots(figsize=(max(10, len(names) * 0.8), 6))
    ax.bar(x - w/2, Si["S1"], w, yerr=Si["S1_conf"], label="S1", capsize=4, color="steelblue")
    ax.bar(x + w/2, Si["ST"], w, yerr=Si["ST_conf"], label="ST", capsize=4, color="firebrick")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=45, ha="right", fontsize=9)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.legend()
    style_axes(ax, title=f"Sobol B — {metric_name}", ylabel="Sobol index")
    plt.tight_layout()
    plt.savefig(filename, dpi=120)
    plt.close()
    print(f"  Saved {filename}")


# -----------------------------------------------------------------------------
# 5. MAIN
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    print("Loading data...")
    df      = load_data("data/raw/paris_weather.csv")
    df_week = df.loc["2024-08-01":"2024-08-07"]   # August — cold circuit fires

    t_array     = np.arange(len(df_week)) * 3600.0
    T_ext_array = df_week["temperature_2m"].values
    dates       = df_week.index
    profiles    = load_profiles("data/raw/Defense_Occupation_Normalised.xlsx")
    Q_int_base, n_ppl_base = build_Q_array(dates, profiles)

    print("Running baseline ODE...")
    T_in_base = run_baseline_ode(t_array, T_ext_array, Q_int_base, n_ppl_base, dates)
    print(f"  T_in: {T_in_base.min():.1f}–{T_in_base.max():.1f}°C")

    N = 512
    samples_B = sobol_sample.sample(PROBLEM_B, N)
    print(f"\nSOBOL B — {PROBLEM_B['num_vars']} params | N={N} | {samples_B.shape[0]} rows")

    N_B = samples_B.shape[0]
    Y_peak_cool  = np.full(N_B, np.nan)
    Y_cumul_cool = np.full(N_B, np.nan)
    Y_peak_heat  = np.full(N_B, np.nan)

    t0 = time.perf_counter()
    for i, row in enumerate(tqdm(samples_B, desc="Sobol B")):
        Y_peak_cool[i], Y_cumul_cool[i], Y_peak_heat[i] = eval_Q_water(
            row, PROBLEM_B["names"], T_in_base, T_ext_array, n_ppl_base, dates)
    print(f"  Done in {time.perf_counter()-t0:.1f}s")

    np.save("data/processed/sobol_B_peak_cool.npy",  Y_peak_cool)
    np.save("data/processed/sobol_B_cumul_cool.npy", Y_cumul_cool)
    np.save("data/processed/sobol_B_peak_heat.npy",  Y_peak_heat)

    Si_pc = sobol_analyze.analyze(PROBLEM_B, Y_peak_cool,  print_to_console=False)
    Si_cc = sobol_analyze.analyze(PROBLEM_B, Y_cumul_cool, print_to_console=False)
    Si_ph = sobol_analyze.analyze(PROBLEM_B, Y_peak_heat,  print_to_console=False)

    plot_sobol(Si_pc, PROBLEM_B["names"], "Peak Q_water_cool [m³/s]",    "images/sobol_B_peak_cool.png")
    plot_sobol(Si_cc, PROBLEM_B["names"], "Cumulative Q_water_cool [m³]", "images/sobol_B_cumul_cool.png")
    plot_sobol(Si_ph, PROBLEM_B["names"], "Peak Q_water_heat [m³/s]",    "images/sobol_B_peak_heat.png")

    print("\nTOP DRIVERS — Peak Q_water_cool (ST)")
    for name, st in sorted(zip(PROBLEM_B["names"], Si_pc["ST"]), key=lambda x: x[1], reverse=True):
        print(f"  {name:<22} ST = {st:.3f}")