# =============================================================================
# sobol.py — Global sensitivity analysis
# 15 operator parameters × 3 output metrics (cost, comfort, CO₂)
# Saltelli sampling (N=512) → 16,384 annual simulations → ~11 hours
# Checkpoint every 500 runs — safe to resume after crash/power loss
# =============================================================================

import os
import time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from SALib.sample import saltelli
from SALib.analyze import sobol

import constants
import regulation
import emissions as em_mod
import humidity as hum_mod
from simulation import run_simulation
from utils import style_axes


# -----------------------------------------------------------------------------
# 1. PROBLEM DEFINITION — 15 operator choices
# -----------------------------------------------------------------------------
# Indices 0–3:  Pareto levers (also swept here for full picture)
# Indices 4–14: remaining operator choices not in Pareto
#
# Ordering constraints:
#   T_DEAD_LOW (4) < T_DEAD_HIGH (5)     — bounds don't overlap, always satisfied
#   T_HW_EXT_HYST (12) < T_HW_EXT_HIGH (11) — bounds overlap, swap if violated
#   T_CW_EXT_HYST (14) > T_CW_EXT_LOW (13)  — bounds overlap, swap if violated

PROBLEM = {
    "num_vars": 15,
    "names": [
        "T_HEAT_LOW_C",           # 0  — heating setpoint [baseline 18]
        "T_COOL_FIXED_C",         # 1  — cooling setpoint [baseline 26]
        "T_HW_SUPPLY_MAX",        # 2  — max hot water supply [baseline 50]
        "AIRFLOW_PER_PERSON_M3H", # 3  — per-person ventilation [baseline 25]
        "T_DEAD_LOW_C",           # 4  — dead band lower T_ext boundary [baseline 15]
        "T_DEAD_HIGH_C",          # 5  — dead band upper T_ext boundary [baseline 22]
        "T_BLOW_HEAT_C",          # 6  — AHU supply air, heating [baseline 30]
        "T_BLOW_COOL_C",          # 7  — AHU supply air, cooling [baseline 15]
        "FRAC_RETURN_AIR",        # 8  — return air fraction [baseline 0.70]
        "T_NIGHT_SETBACK_C",      # 9  — anti-freeze setpoint [baseline 5]
        "T_STAIR_COLD_C",         # 10 — curtain activation T_ext [baseline 7]
        "T_HW_EXT_HIGH_C",        # 11 — hot water shutoff T_ext [baseline 15]
        "T_HW_EXT_HYST_C",        # 12 — hot water restart T_ext [baseline 13]
        "T_CW_EXT_LOW_C",         # 13 — cold water shutoff T_ext [baseline 26]
        "T_CW_EXT_HYST_C",        # 14 — cold water restart T_ext [baseline 27]
    ],
    "bounds": [
        [16, 22],       # 0  T_HEAT_LOW_C
        [24, 28],       # 1  T_COOL_FIXED_C
        [40, 55],       # 2  T_HW_SUPPLY_MAX
        [15, 35],       # 3  AIRFLOW_PER_PERSON_M3H
        [12, 18],       # 4  T_DEAD_LOW_C
        [20, 25],       # 5  T_DEAD_HIGH_C        (always > idx 4)
        [28, 35],       # 6  T_BLOW_HEAT_C
        [13, 17],       # 7  T_BLOW_COOL_C
        [0.50, 0.85],   # 8  FRAC_RETURN_AIR
        [3, 8],         # 9  T_NIGHT_SETBACK_C
        [3, 10],        # 10 T_STAIR_COLD_C
        [13, 18],       # 11 T_HW_EXT_HIGH_C
        [10, 15],       # 12 T_HW_EXT_HYST_C      (overlap with 11: swap if > 11)
        [24, 28],       # 13 T_CW_EXT_LOW_C
        [25, 29],       # 14 T_CW_EXT_HYST_C      (overlap with 13: swap if < 13)
    ],
}

N_SAMPLES = 512
N_TOTAL = N_SAMPLES * (2 * PROBLEM["num_vars"] + 2)   # 16,384

CHECKPOINT_PATH = os.path.join("data", "processed", "sobol_checkpoint.npz")
RESULTS_PATH    = os.path.join("data", "processed", "sobol_results.npz")


# -----------------------------------------------------------------------------
# 2. MONKEY-PATCH HELPERS
# -----------------------------------------------------------------------------

_ORIGINALS = {}


def _patch(name: str, value: float):
    """Patch a constant in constants + every module that imported it by name."""
    if name not in _ORIGINALS:
        _ORIGINALS[name] = getattr(constants, name)

    setattr(constants, name, value)
    for mod in (regulation, em_mod, hum_mod):
        if hasattr(mod, name):
            setattr(mod, name, value)


def _patch_derived():
    """Recompute constants that depend on swept values."""
    ov = constants.AIRFLOW_OVERPRESSURE_M3H
    pp = constants.AIRFLOW_PER_PERSON_M3H
    pk = constants.PEOPLE_PEAK

    _patch("AIRFLOW_MIN_M3H", ov)
    _patch("AIRFLOW_MAX_M3H", (ov + pk * pp) * 1.10)
    _patch("P_FAN_RATED_W",
           constants.AIRFLOW_MAX_M3H / 3600.0 * constants.DP_AHU_PA / constants.ETA_FAN)

    _patch("T_HEAT_HIGH_C", constants.T_HEAT_LOW_C + 2.0)


def _restore_all():
    """Restore all constants to original values."""
    for name, val in _ORIGINALS.items():
        setattr(constants, name, val)
        for mod in (regulation, em_mod, hum_mod):
            if hasattr(mod, name):
                setattr(mod, name, val)


# -----------------------------------------------------------------------------
# 3. ORDERING CONSTRAINTS
# -----------------------------------------------------------------------------

def _fix_ordering(row: np.ndarray) -> np.ndarray:
    """Swap values when Saltelli sampling violates physical ordering."""
    # T_HW_EXT_HYST (12) must be < T_HW_EXT_HIGH (11)
    if row[12] > row[11]:
        row[11], row[12] = row[12], row[11]

    # T_CW_EXT_HYST (14) must be > T_CW_EXT_LOW (13)
    if row[14] < row[13]:
        row[13], row[14] = row[14], row[13]

    return row


# -----------------------------------------------------------------------------
# 4. SINGLE MODEL EVALUATION
# -----------------------------------------------------------------------------

def eval_row(row: np.ndarray) -> tuple:
    """
    Patch all 15 params from one Saltelli row, run annual sim,
    return (cost_eur, combined_discomfort_pct, CO2_kgCO2).
    """
    row = _fix_ordering(row.copy())

    names = PROBLEM["names"]
    for j, name in enumerate(names):
        _patch(name, float(row[j]))
    _patch_derived()

    try:
        r = run_simulation()
    except Exception as e:
        print(f"  FAILED: {e}")
        return (np.nan, np.nan, np.nan)

    cost    = r["em"]["cost_annual_eur"]
    comfort = r["comfort"]["combined_discomfort_pct"]
    co2     = r["em"]["CO2_annual_kgCO2"]

    return (cost, comfort, co2)


# -----------------------------------------------------------------------------
# 5. PLOTTING
# -----------------------------------------------------------------------------

def plot_sobol(Si: dict, metric_name: str, metric_unit: str,
               names: list, filename: str):
    """Grouped bar chart: S1 and ST with confidence intervals, sorted by ST."""
    s1 = np.array(Si["S1"])
    st = np.array(Si["ST"])
    s1_conf = np.array(Si["S1_conf"])
    st_conf = np.array(Si["ST_conf"])

    # Sort by ST descending
    order = np.argsort(st)[::-1]
    s1       = s1[order]
    st       = st[order]
    s1_conf  = s1_conf[order]
    st_conf  = st_conf[order]
    sorted_names = [names[i] for i in order]

    n = len(names)
    x = np.arange(n)
    w = 0.35

    fig, ax = plt.subplots(figsize=(14, 6))

    ax.bar(x - w/2, s1, w, yerr=s1_conf, capsize=3,
           label="S1 (first-order)", color="#3498db", alpha=0.85)
    ax.bar(x + w/2, st, w, yerr=st_conf, capsize=3,
           label="ST (total-order)", color="#e67e22", alpha=0.85)

    ax.set_xticks(x)
    ax.set_xticklabels(sorted_names, rotation=45, ha="right", fontsize=9)
    ax.axhline(0, color="gray", lw=0.5)

    style_axes(ax,
               title=f"Sobol Sensitivity — {metric_name} [{metric_unit}]",
               ylabel="Sensitivity index")
    ax.legend(fontsize=10, loc="upper right")

    sum_s1 = Si["S1"].sum()
    ax.text(0.99, 0.85, f"Σ S1 = {sum_s1:.2f}",
            transform=ax.transAxes, ha="right", fontsize=10,
            bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Saved: {filename}")


# -----------------------------------------------------------------------------
# 6. RESULTS TABLE
# -----------------------------------------------------------------------------

def print_results(Si: dict, metric_name: str, names: list):
    """Print S1 / ST table sorted by ST, with confidence intervals."""
    s1   = Si["S1"]
    st   = Si["ST"]
    s1_c = Si["S1_conf"]
    st_c = Si["ST_conf"]

    order = np.argsort(st)[::-1]

    print(f"\n{'='*60}")
    print(f" {metric_name}")
    print(f"{'='*60}")
    print(f"{'Parameter':<28} {'S1':>7} {'±':>5} {'ST':>7} {'±':>5}")
    print(f"{'-'*60}")
    for i in order:
        print(f"  {names[i]:<26} {s1[i]:>7.3f} {s1_c[i]:>5.3f} {st[i]:>7.3f} {st_c[i]:>5.3f}")
    print(f"{'-'*60}")
    print(f"  {'Σ S1':<26} {s1.sum():>7.3f}")
    print()


# -----------------------------------------------------------------------------
# 7. MAIN
# -----------------------------------------------------------------------------

if __name__ == "__main__":

    print(f"Sobol GSA: {PROBLEM['num_vars']} parameters, N={N_SAMPLES}")
    print(f"Total runs: {N_TOTAL}")
    print(f"Estimated time: {N_TOTAL * 2.4 / 3600:.1f} hours\n")

    # --- Generate Saltelli samples ---
    X = saltelli.sample(PROBLEM, N_SAMPLES)
    assert X.shape == (N_TOTAL, PROBLEM["num_vars"])

    # --- Allocate output arrays ---
    Y_cost    = np.full(N_TOTAL, np.nan)
    Y_comfort = np.full(N_TOTAL, np.nan)
    Y_co2     = np.full(N_TOTAL, np.nan)

    # --- Resume from checkpoint if available ---
    start = 0
    if os.path.exists(CHECKPOINT_PATH):
        ckpt = np.load(CHECKPOINT_PATH)
        start = int(ckpt["completed"])
        Y_cost[:start]    = ckpt["Y_cost"]
        Y_comfort[:start] = ckpt["Y_comfort"]
        Y_co2[:start]     = ckpt["Y_co2"]
        print(f"Resumed from checkpoint: {start}/{N_TOTAL} completed\n")

    # --- Run loop ---
    t0 = time.perf_counter()

    for i in range(start, N_TOTAL):
        cost, comfort, co2 = eval_row(X[i])
        Y_cost[i]    = cost
        Y_comfort[i] = comfort
        Y_co2[i]     = co2

        # Progress every 100 runs
        if (i + 1) % 100 == 0 or i == start:
            elapsed = time.perf_counter() - t0
            done = i + 1 - start
            rate = done / elapsed if elapsed > 0 else 0
            remaining = (N_TOTAL - i - 1) / rate / 3600 if rate > 0 else 0
            print(f"  [{i+1}/{N_TOTAL}]  {rate:.1f} runs/s  "
                  f"ETA {remaining:.1f}h  "
                  f"cost={cost:.0f}€  comfort={comfort:.1f}%  co2={co2:.0f}kg")

        # Checkpoint every 500 runs
        if (i + 1) % 500 == 0:
            np.savez(CHECKPOINT_PATH,
                     Y_cost=Y_cost[:i+1],
                     Y_comfort=Y_comfort[:i+1],
                     Y_co2=Y_co2[:i+1],
                     completed=i+1)

    _restore_all()

    elapsed_total = time.perf_counter() - t0
    print(f"\nCompleted {N_TOTAL} runs in {elapsed_total/3600:.1f}h "
          f"({N_TOTAL/elapsed_total:.1f} runs/s)")

    # --- Save full results ---
    np.savez(RESULTS_PATH, X=X, Y_cost=Y_cost, Y_comfort=Y_comfort, Y_co2=Y_co2)
    print(f"Saved: {RESULTS_PATH}")

    # --- Check for failed runs ---
    n_nan = np.isnan(Y_cost).sum()
    if n_nan > 0:
        print(f"WARNING: {n_nan} runs failed (NaN). Sobol indices may be unreliable.")

    # --- Analyze ---
    names = PROBLEM["names"]

    Si_cost    = sobol.analyze(PROBLEM, Y_cost,    print_to_console=False)
    Si_comfort = sobol.analyze(PROBLEM, Y_comfort, print_to_console=False)
    Si_co2     = sobol.analyze(PROBLEM, Y_co2,     print_to_console=False)

    # --- Print tables ---
    print_results(Si_cost,    "Annual Cost (€)",          names)
    print_results(Si_comfort, "Combined Discomfort (%)",   names)
    print_results(Si_co2,     "Annual CO₂ (kgCO₂)",       names)

    # --- Plot ---
    plot_sobol(Si_cost,    "Annual Cost",          "€",     names, "images/sobol_cost.png")
    plot_sobol(Si_comfort, "Combined Discomfort",  "%",     names, "images/sobol_comfort.png")
    plot_sobol(Si_co2,     "Annual CO₂",           "kgCO₂", names, "images/sobol_co2.png")

    # --- Cleanup checkpoint ---
    if os.path.exists(CHECKPOINT_PATH):
        os.remove(CHECKPOINT_PATH)
        print("Checkpoint removed (run complete).")
