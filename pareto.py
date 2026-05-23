# =============================================================================
# pareto.py — Multi-criteria Pareto sweep over regulation setpoints
# Axes: E_total (kWh/year) vs discomfort (% service hours outside 18–26°C)
# Control levers: T_HEAT_LOW, T_COOL_FIXED, T_HW_SUPPLY_MAX,
#                 AIRFLOW_OVERPRESSURE, AIRFLOW_PER_PERSON
# =============================================================================

import time
import itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import constants
import regulation
import emissions as em_mod
from thermal_model import run_simulation


# -----------------------------------------------------------------------------
# 1. PARAMETER GRID
# -----------------------------------------------------------------------------

LEVERS = {
    "T_HEAT_LOW_C":            np.array([16.0, 18.0, 20.0, 21.0]),
    "T_COOL_FIXED_C":          np.array([24.0, 25.0, 26.0, 28.0]),
    "T_HW_SUPPLY_MAX":         np.array([40.0, 45.0, 50.0, 55.0]),
    "AIRFLOW_OVERPRESSURE_M3H": np.array([1500.0, 2500.0, 3500.0]),
    "AIRFLOW_PER_PERSON_M3H":  np.array([15.0, 25.0, 35.0]),
}

# Total configs: 4 × 4 × 4 × 3 × 3 = 576
N_CONFIGS = 1
for v in LEVERS.values():
    N_CONFIGS *= len(v)


# -----------------------------------------------------------------------------
# 2. MONKEY-PATCH HELPER
# -----------------------------------------------------------------------------

# Save originals to restore at the end
_ORIGINALS = {}

def _patch(name: str, value: float):
    """Patch a constant in constants.py AND in every module that imported it by name."""
    if name not in _ORIGINALS:
        _ORIGINALS[name] = getattr(constants, name)

    setattr(constants, name, value)
    # regulation.py and emissions.py use `from constants import X`
    # so the local name must be patched too
    if hasattr(regulation, name):
        setattr(regulation, name, value)
    if hasattr(em_mod, name):
        setattr(em_mod, name, value)


def _patch_derived():
    """Recompute derived constants that depend on swept values."""
    ov = constants.AIRFLOW_OVERPRESSURE_M3H
    pp = constants.AIRFLOW_PER_PERSON_M3H
    pk = constants.PEOPLE_PEAK

    new_min = ov
    new_max = (ov + pk * pp) * 1.10
    new_pfan = new_max / 3600.0 * constants.DP_AHU_PA / constants.ETA_FAN

    _patch("AIRFLOW_MIN_M3H", new_min)
    _patch("AIRFLOW_MAX_M3H", new_max)
    _patch("P_FAN_RATED_W", new_pfan)

    # T_HEAT_HIGH follows T_HEAT_LOW: ramp ends 2°C above
    _patch("T_HEAT_HIGH_C", constants.T_HEAT_LOW_C + 2.0)


def _restore_all():
    """Restore all constants to original values."""
    for name, val in _ORIGINALS.items():
        setattr(constants, name, val)
        if hasattr(regulation, name):
            setattr(regulation, name, val)
        if hasattr(em_mod, name):
            setattr(em_mod, name, val)


# -----------------------------------------------------------------------------
# 3. SINGLE EVALUATION
# -----------------------------------------------------------------------------

def eval_config(config: dict) -> dict | None:
    """
    Patch constants, run annual sim, extract objectives.
    Returns dict with config + objectives, or None on failure.
    """
    # Patch all levers
    for name, val in config.items():
        _patch(name, val)
    _patch_derived()

    try:
        r = run_simulation()
    except Exception as e:
        print(f"  FAILED: {e}")
        return None

    em = r["em"]
    T_in = r["T_in"]
    dates = r["dates"]

    # Service hours mask
    is_service = np.array([not (1 <= ts.hour < 5) for ts in dates])
    n_service = is_service.sum()

    pct_above_26 = ((T_in > 26.0) & is_service).sum() / n_service * 100.0
    pct_below_18 = ((T_in < 18.0) & is_service).sum() / n_service * 100.0
    discomfort = pct_above_26 + pct_below_18

    return {
        **config,
        "E_total_kWh": em["E_annual_kWh"],
        "E_heat_kWh":  em["E_heat_total_kWh"],
        "E_cool_kWh":  em["E_cool_total_kWh"],
        "E_fan_kWh":   em["E_fan_total_kWh"],
        "E_curtain_kWh": em["E_curtain_total_kWh"],
        "cost_eur":    em["cost_annual_eur"],
        "CO2_kg":      em["CO2_annual_kgCO2"],
        "pct_above_26": pct_above_26,
        "pct_below_18": pct_below_18,
        "discomfort":  discomfort,
        "T_in_min":    T_in.min(),
        "T_in_max":    T_in.max(),
    }


# -----------------------------------------------------------------------------
# 4. PARETO FRONT EXTRACTION
# -----------------------------------------------------------------------------

def pareto_front(points: list[dict], x_key="E_total_kWh", y_key="discomfort") -> list[int]:
    """Return indices of non-dominated points (minimize both axes)."""
    n = len(points)
    is_dominated = [False] * n
    for i in range(n):
        if is_dominated[i]:
            continue
        for j in range(n):
            if i == j or is_dominated[j]:
                continue
            # j dominates i if j ≤ i on both and j < i on at least one
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
    """Scatter all configs, highlight Pareto front."""
    E_all = [r["E_total_kWh"] for r in results]
    D_all = [r["discomfort"] for r in results]

    E_front = [results[i]["E_total_kWh"] for i in front_idx]
    D_front = [results[i]["discomfort"] for i in front_idx]

    # Sort front by energy for the line
    order = np.argsort(E_front)
    E_front_sorted = [E_front[i] for i in order]
    D_front_sorted = [D_front[i] for i in order]

    fig, ax = plt.subplots(figsize=(12, 8))

    ax.scatter(E_all, D_all, s=18, alpha=0.3, color="gray", label="Dominated")
    ax.scatter(E_front_sorted, D_front_sorted, s=60, color="crimson",
               zorder=5, label="Pareto front")
    ax.plot(E_front_sorted, D_front_sorted, color="crimson", lw=1.5,
            ls="--", alpha=0.7, zorder=4)

    # Annotate a few front points
    for idx in front_idx:
        r = results[idx]
        label = (f"H={r['T_HEAT_LOW_C']:.0f} C={r['T_COOL_FIXED_C']:.0f} "
                 f"HW={r['T_HW_SUPPLY_MAX']:.0f}")
        ax.annotate(label, (r["E_total_kWh"], r["discomfort"]),
                    fontsize=6, alpha=0.7,
                    xytext=(5, 5), textcoords="offset points")

    ax.set_xlabel("Annual energy [kWh]", fontsize=12)
    ax.set_ylabel("Discomfort [% service hours outside 18–26°C]", fontsize=12)
    ax.set_title(f"Pareto front — {len(results)} configurations, "
                 f"{len(front_idx)} non-dominated",
                 fontsize=13)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Saved: {filename}")


# -----------------------------------------------------------------------------
# 6. MAIN
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    # Build all combinations
    lever_names = list(LEVERS.keys())
    lever_values = list(LEVERS.values())
    all_combos = list(itertools.product(*lever_values))

    print(f"Pareto sweep: {N_CONFIGS} configurations, {len(lever_names)} levers")
    print(f"Levers: {lever_names}")

    results = []
    t0 = time.perf_counter()

    for i, combo in enumerate(all_combos):
        config = dict(zip(lever_names, combo))
        if (i + 1) % 50 == 0 or i == 0:
            elapsed = time.perf_counter() - t0
            rate = (i + 1) / elapsed if elapsed > 0 else 0
            eta = (N_CONFIGS - i - 1) / rate / 60 if rate > 0 else 0
            print(f"  [{i+1}/{N_CONFIGS}] {rate:.1f} runs/min, ETA {eta:.0f} min — {config}")

        result = eval_config(config)
        if result is not None:
            results.append(result)

    _restore_all()

    elapsed = time.perf_counter() - t0
    print(f"\nCompleted {len(results)}/{N_CONFIGS} in {elapsed/60:.1f} min "
          f"({len(results)/elapsed*60:.0f} runs/min)")

    # Find Pareto front
    front_idx = pareto_front(results)
    print(f"Pareto front: {len(front_idx)} non-dominated points")

    # Sort front by energy
    front_idx_sorted = sorted(front_idx, key=lambda i: results[i]["E_total_kWh"])

    print("\n=== Pareto front (sorted by energy) ===")
    print(f"{'E_total':>8} {'discomf':>8} {'%>26':>6} {'%<18':>6} "
          f"{'T_heat':>6} {'T_cool':>6} {'HW_max':>6} {'Ovpr':>6} {'A/pp':>6} {'cost':>6}")
    for idx in front_idx_sorted:
        r = results[idx]
        print(f"{r['E_total_kWh']:>8.0f} {r['discomfort']:>8.1f} "
              f"{r['pct_above_26']:>6.1f} {r['pct_below_18']:>6.1f} "
              f"{r['T_HEAT_LOW_C']:>6.0f} {r['T_COOL_FIXED_C']:>6.0f} "
              f"{r['T_HW_SUPPLY_MAX']:>6.0f} "
              f"{r['AIRFLOW_OVERPRESSURE_M3H']:>6.0f} "
              f"{r['AIRFLOW_PER_PERSON_M3H']:>6.0f} "
              f"{r['cost_eur']:>6.0f}")

    plot_pareto(results, front_idx)