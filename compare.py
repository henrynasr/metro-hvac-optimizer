# =============================================================================
# compare.py — Baseline vs Pareto-optimal configuration
# Runs both configs, prints side-by-side comparison, saves plots for each.
# Justification: Sobol GSA (15 params) → Pareto front (1,296 configs, 6 levers)
# =============================================================================

import numpy as np
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
# 1. OPTIMAL CONFIGURATION — from Pareto front (cheapest at H=16)
# -----------------------------------------------------------------------------

OPTIMAL = {
    "T_HEAT_LOW_C":           16.0,    # was 18.0 — Sobol S1=0.86 on cost
    "AIRFLOW_PER_PERSON_M3H": 18.0,    # was 25.0 — minimum regulatory
    "T_HW_SUPPLY_MAX":        40.0,    # was 50.0 — lower = better COP, always
    "T_STAIR_COLD_C":          7.0,    # unchanged — curtain at 7°C
    "T_HW_EXT_HYST_C":       13.0,    # unchanged — no effect on front
    "T_BLOW_HEAT_C":          28.0,    # was 30.0 — cheaper, similar comfort
}


# -----------------------------------------------------------------------------
# 2. MONKEY-PATCH HELPERS
# -----------------------------------------------------------------------------

_ORIGINALS = {}


def _patch(name: str, value: float):
    if name not in _ORIGINALS:
        _ORIGINALS[name] = getattr(constants, name)
    setattr(constants, name, value)
    for mod in (regulation, em_mod, hum_mod):
        if hasattr(mod, name):
            setattr(mod, name, value)


def _patch_derived():
    ov = constants.AIRFLOW_OVERPRESSURE_M3H
    pp = constants.AIRFLOW_PER_PERSON_M3H
    pk = constants.PEOPLE_PEAK
    _patch("AIRFLOW_MIN_M3H", ov)
    _patch("AIRFLOW_MAX_M3H", (ov + pk * pp) * 1.10)
    _patch("P_FAN_RATED_W",
           constants.AIRFLOW_MAX_M3H / 3600.0 * constants.DP_AHU_PA / constants.ETA_FAN)
    _patch("T_HEAT_HIGH_C", constants.T_HEAT_LOW_C + 2.0)


def _restore_all():
    for name, val in _ORIGINALS.items():
        setattr(constants, name, val)
        for mod in (regulation, em_mod, hum_mod):
            if hasattr(mod, name):
                setattr(mod, name, val)


# -----------------------------------------------------------------------------
# 3. RUN AND EXTRACT
# -----------------------------------------------------------------------------

def extract_summary(r: dict) -> dict:
    """Extract key metrics from simulation results."""
    em = r["em"]
    c  = r["comfort"]
    return {
        "E_total_kWh":     em["E_annual_kWh"],
        "E_heat_kWh":      em["E_heat_total_kWh"],
        "E_cool_kWh":      em["E_cool_total_kWh"],
        "E_fan_kWh":       em["E_fan_total_kWh"],
        "E_curtain_kWh":   em["E_curtain_total_kWh"],
        "cost_eur":        em["cost_annual_eur"],
        "CO2_kg":          em["CO2_annual_kgCO2"],
        "T_in_min":        r["T_in"].min(),
        "T_in_max":        r["T_in"].max(),
        "T_comfort":       c["T_comfort_pct"],
        "T_mild":          c["T_mild_pct"],
        "T_discomfort":    c["T_discomfort_pct"],
        "RH_comfort":      c["RH_comfort_pct"],
        "RH_mild":         c["RH_mild_pct"],
        "RH_discomfort":   c["RH_discomfort_pct"],
        "combined_comfort":    c["combined_comfort_pct"],
        "combined_mild":       c["combined_mild_pct"],
        "combined_discomfort": c["combined_discomfort_pct"],
        "condensation_h":  r["humidity"]["hours_condensation"],
        "curtain_h":       int(r["curtain_on"].sum()),
    }


# -----------------------------------------------------------------------------
# 4. PRINT COMPARISON
# -----------------------------------------------------------------------------

def print_comparison(base: dict, opt: dict):
    """Side-by-side table with deltas."""

    def row(label, key, unit="", fmt=".0f", invert=False):
        b = base[key]
        o = opt[key]
        d = o - b
        sign = "+" if d > 0 else ""
        if abs(d) < 0.05:
            delta = "—"
        elif unit == "%":
            delta = f"{sign}{d:{fmt}}pp"
        else:
            pct = d / b * 100 if b != 0 else 0
            delta = f"{sign}{d:{fmt}}{unit} ({sign}{pct:.0f}%)"
        print(f"  {label:<28} {b:>10{fmt}}{unit}  {o:>10{fmt}}{unit}  {delta}")

    print()
    print(f"{'='*75}")
    print(f" Baseline vs Pareto-Optimal — Annual 2024")
    print(f"{'='*75}")

    # Config differences
    print(f"\n  {'Parameter':<28} {'Baseline':>10}  {'Optimal':>10}  {'Change'}")
    print(f"  {'-'*70}")
    for name, opt_val in OPTIMAL.items():
        base_val = _ORIGINALS.get(name, getattr(constants, name))
        if abs(opt_val - base_val) < 0.01:
            delta = "—"
        else:
            delta = f"{opt_val - base_val:+.0f}"
        print(f"  {name:<28} {base_val:>10.0f}  {opt_val:>10.0f}  {delta}")

    # Energy
    print(f"\n  {'Metric':<28} {'Baseline':>10}  {'Optimal':>10}  {'Delta'}")
    print(f"  {'-'*70}")
    row("Heating",             "E_heat_kWh",    " kWh")
    row("Cooling",             "E_cool_kWh",    " kWh")
    row("Fans",                "E_fan_kWh",     " kWh")
    row("Curtain",             "E_curtain_kWh", " kWh")
    row("TOTAL energy",        "E_total_kWh",   " kWh")
    row("TOTAL cost",          "cost_eur",      " €")
    row("TOTAL CO₂",           "CO2_kg",        " kg")

    # Temperature
    print(f"\n  {'Temperature comfort':<28} {'Baseline':>10}  {'Optimal':>10}  {'Delta'}")
    print(f"  {'-'*70}")
    row("T_in range min",      "T_in_min",      "°C", ".1f")
    row("T_in range max",      "T_in_max",      "°C", ".1f")
    row("Comfort (18–26°C)",   "T_comfort",     "%", ".1f")
    row("Mild (14–18/26–28°C)","T_mild",        "%", ".1f")
    row("Discomfort (<14/>28)","T_discomfort",   "%", ".1f")

    # Humidity
    print(f"\n  {'Humidity comfort':<28} {'Baseline':>10}  {'Optimal':>10}  {'Delta'}")
    print(f"  {'-'*70}")
    row("Comfort (40–60%)",    "RH_comfort",    "%", ".1f")
    row("Mild (35–40/60–65%)", "RH_mild",       "%", ".1f")
    row("Discomfort (<35/>65%)","RH_discomfort", "%", ".1f")
    row("Condensation hours",  "condensation_h", " h")

    # Combined
    print(f"\n  {'Combined comfort':<28} {'Baseline':>10}  {'Optimal':>10}  {'Delta'}")
    print(f"  {'-'*70}")
    row("Comfort (both OK)",   "combined_comfort",    "%", ".1f")
    row("Mild (tolerable)",    "combined_mild",       "%", ".1f")
    row("Discomfort (≥1 bad)", "combined_discomfort",  "%", ".1f")
    row("Curtain runtime",     "curtain_h",           " h")

    print()


# -----------------------------------------------------------------------------
# 5. MAIN
# -----------------------------------------------------------------------------

if __name__ == "__main__":

    # --- Run baseline ---
    print("Running baseline...")
    r_base = run_simulation()
    s_base = extract_summary(r_base)

    # --- Run optimal ---
    print("Running optimal...")
    for name, val in OPTIMAL.items():
        _patch(name, val)
    _patch_derived()

    r_opt = run_simulation()
    s_opt = extract_summary(r_opt)

    _restore_all()

    # --- Print comparison ---
    print_comparison(s_base, s_opt)
