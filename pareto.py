# =============================================================================
# pareto.py — Multi-criteria Pareto sweep
# Axes: cost (€/year) vs combined discomfort (% service hours)
# 6 levers: T_HEAT_LOW, AIRFLOW/PP, T_HW_SUPPLY_MAX,
#           T_STAIR_COLD, T_HW_EXT_HYST, T_BLOW_HEAT
# =============================================================================

import time
import itertools
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import constants
import regulation
import emissions as em_mod
import humidity as hum_mod
from simulation import run_simulation
from utils import style_axes


# -----------------------------------------------------------------------------
# 1. PARAMETER GRID — 3×3×4×4×3×3 = 1,296 configs
# -----------------------------------------------------------------------------

LEVERS = {
    "T_HEAT_LOW_C":             np.array([14.0, 16.0, 18.0]),
    "AIRFLOW_PER_PERSON_M3H":   np.array([18.0, 25.0, 30.0]),
    "T_HW_SUPPLY_MAX":          np.array([40.0, 45.0, 50.0, 55.0]),
    "T_STAIR_COLD_C":           np.array([3.0, 5.0, 7.0, 10.0]),
    "T_HW_EXT_HYST_C":          np.array([11.0, 13.0, 15.0]),
    "T_BLOW_HEAT_C":             np.array([28.0, 30.0, 32.0]),
}

N_CONFIGS = 1
for v in LEVERS.values():
    N_CONFIGS *= len(v)


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
# 3. SINGLE EVALUATION
# -----------------------------------------------------------------------------

def eval_config(config: dict) -> dict | None:
    """Patch constants, run annual sim, extract objectives."""
    for name, val in config.items():
        _patch(name, val)
    _patch_derived()

    try:
        r = run_simulation()
    except Exception as e:
        print(f"  FAILED: {e}")
        return None

    em = r["em"]
    c  = r["comfort"]

    return {
        **config,
        # Energy
        "E_total_kWh":    em["E_annual_kWh"],
        "E_heat_kWh":     em["E_heat_total_kWh"],
        "E_cool_kWh":     em["E_cool_total_kWh"],
        "E_fan_kWh":      em["E_fan_total_kWh"],
        "E_curtain_kWh":  em["E_curtain_total_kWh"],
        "cost_eur":       em["cost_annual_eur"],
        "CO2_kg":         em["CO2_annual_kgCO2"],
        # Temperature comfort
        "T_comfort_pct":    c["T_comfort_pct"],
        "T_mild_pct":       c["T_mild_pct"],
        "T_discomfort_pct": c["T_discomfort_pct"],
        # Humidity comfort
        "RH_comfort_pct":    c["RH_comfort_pct"],
        "RH_mild_pct":       c["RH_mild_pct"],
        "RH_discomfort_pct": c["RH_discomfort_pct"],
        # Combined comfort
        "combined_comfort_pct":    c["combined_comfort_pct"],
        "combined_mild_pct":       c["combined_mild_pct"],
        "combined_discomfort_pct": c["combined_discomfort_pct"],
        # Ranges
        "T_in_min": r["T_in"].min(),
        "T_in_max": r["T_in"].max(),
    }


# -----------------------------------------------------------------------------
# 4. PARETO FRONT EXTRACTION
# -----------------------------------------------------------------------------

def pareto_front(points: list[dict],
                 x_key="cost_eur",
                 y_key="combined_discomfort_pct") -> list[int]:
    """Return indices of non-dominated points (minimize both axes)."""
    n = len(points)
    is_dominated = [False] * n
    for i in range(n):
        if is_dominated[i]:
            continue
        for j in range(n):
            if i == j or is_dominated[j]:
                continue
            if (points[j][x_key] <= points[i][x_key] and
                points[j][y_key] <= points[i][y_key] and
                (points[j][x_key] < points[i][x_key] or
                 points[j][y_key] < points[i][y_key])):
                is_dominated[i] = True
                break
    return [i for i in range(n) if not is_dominated[i]]


# -----------------------------------------------------------------------------
# 5. PLOT
# -----------------------------------------------------------------------------

def plot_pareto(results: list[dict], front_idx: list[int],
                filename: str = "images/pareto_front.png"):
    """Scatter all configs, highlight Pareto front, label key points only."""
    x_all = [r["cost_eur"] for r in results]
    y_all = [r["combined_discomfort_pct"] for r in results]

    x_front = [results[i]["cost_eur"] for i in front_idx]
    y_front = [results[i]["combined_discomfort_pct"] for i in front_idx]

    # Sort front by cost for the connecting line
    order = np.argsort(x_front)
    x_sorted = [x_front[i] for i in order]
    y_sorted = [y_front[i] for i in order]
    idx_sorted = [front_idx[i] for i in order]

    fig, ax = plt.subplots(figsize=(14, 8))

    ax.scatter(x_all, y_all, s=12, alpha=0.2, color="gray", label="Dominated")
    ax.scatter(x_sorted, y_sorted, s=50, color="crimson",
               zorder=5, label="Pareto front")
    ax.plot(x_sorted, y_sorted, color="crimson", lw=1.5,
            ls="--", alpha=0.6, zorder=4)

    # Label only every 4th front point to avoid clutter
    for k, idx in enumerate(idx_sorted):
        if k % 4 != 0 and k != len(idx_sorted) - 1:
            continue
        r = results[idx]
        label = (f"H={r['T_HEAT_LOW_C']:.0f} "
                 f"A={r['AIRFLOW_PER_PERSON_M3H']:.0f} "
                 f"HW={r['T_HW_SUPPLY_MAX']:.0f}")
        ax.annotate(label, (r["cost_eur"], r["combined_discomfort_pct"]),
                    fontsize=7, alpha=0.8,
                    xytext=(8, 5), textcoords="offset points",
                    bbox=dict(boxstyle="round,pad=0.2",
                              facecolor="white", alpha=0.7, edgecolor="none"))

    style_axes(ax,
               title=f"Pareto Front — {len(results)} configs, "
                     f"{len(front_idx)} non-dominated",
               xlabel="Annual cost [€]",
               ylabel="Combined discomfort [% service hours]")
    ax.legend(fontsize=10, loc="upper right")

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Saved: {filename}")


# -----------------------------------------------------------------------------
# 6. MAIN
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    lever_names = list(LEVERS.keys())
    lever_values = list(LEVERS.values())
    all_combos = list(itertools.product(*lever_values))

    print(f"Pareto sweep: {N_CONFIGS} configurations, {len(lever_names)} levers")
    for name, vals in LEVERS.items():
        print(f"  {name}: {vals}")

    results = []
    t0 = time.perf_counter()

    for i, combo in enumerate(all_combos):
        config = dict(zip(lever_names, combo))
        if (i + 1) % 100 == 0 or i == 0:
            elapsed = time.perf_counter() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (N_CONFIGS - i - 1) / rate / 60 if rate > 0 else 0
            print(f"  [{i+1}/{N_CONFIGS}] {rate:.1f} runs/min "
                  f"ETA {eta:.0f} min — H={config['T_HEAT_LOW_C']:.0f} "
                  f"A={config['AIRFLOW_PER_PERSON_M3H']:.0f} "
                  f"HW={config['T_HW_SUPPLY_MAX']:.0f}")

        result = eval_config(config)
        if result is not None:
            results.append(result)

    _restore_all()
    
    # --- Save all results ---
    df_results = pd.DataFrame(results)
    df_results.to_csv("data/processed/pareto_all_configs.csv", index=False)
    print(f"Saved: data/processed/pareto_all_configs.csv ({len(df_results)} rows)")

    elapsed = time.perf_counter() - t0
    print(f"\nCompleted {len(results)}/{N_CONFIGS} in {elapsed/60:.1f} min "
          f"({len(results)/elapsed*60:.0f} runs/min)")

    # --- Pareto front ---
    front_idx = pareto_front(results)
    front_idx_sorted = sorted(front_idx, key=lambda i: results[i]["cost_eur"])

    print(f"Pareto front: {len(front_idx)} non-dominated points\n")

    # --- Print table ---
    print(f"{'cost':>7} {'discomf':>7} {'T_comf':>6} {'T_mild':>6} {'T_disc':>6} "
          f"{'RH_comf':>7} {'RH_dis':>6} "
          f"{'H':>4} {'A/pp':>4} {'HW':>4} {'Stair':>5} {'Hyst':>4} {'Blow':>4}")
    print("-" * 95)

    for idx in front_idx_sorted:
        r = results[idx]
        print(f"{r['cost_eur']:>7.0f} {r['combined_discomfort_pct']:>7.1f} "
              f"{r['T_comfort_pct']:>6.1f} {r['T_mild_pct']:>6.1f} "
              f"{r['T_discomfort_pct']:>6.1f} "
              f"{r['RH_comfort_pct']:>7.1f} {r['RH_discomfort_pct']:>6.1f} "
              f"{r['T_HEAT_LOW_C']:>4.0f} "
              f"{r['AIRFLOW_PER_PERSON_M3H']:>4.0f} "
              f"{r['T_HW_SUPPLY_MAX']:>4.0f} "
              f"{r['T_STAIR_COLD_C']:>5.0f} "
              f"{r['T_HW_EXT_HYST_C']:>4.0f} "
              f"{r['T_BLOW_HEAT_C']:>4.0f}")

    # --- Plot ---
    plot_pareto(results, front_idx)