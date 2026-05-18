# =============================================================================
# occupancy.py — Occupancy profiles and infiltration rate, one platform zone
# Functions: load_profiles, build_Q_array, v_inf_m3s, _day_type
# =============================================================================

import numpy as np
import pandas as pd

from constants import (
    PEOPLE_PEAK, WATTS_SENSIBLE_PP, BASELINE_W,
    ETA_INF, V_CYCLE_M3,
    HEADWAY_PEAK_S, HEADWAY_OFFPEAK_S, HEADWAY_NIGHT_S,
)


# -----------------------------------------------------------------------------
# 1. HOLIDAY CALENDAR — Zone C 2024
# -----------------------------------------------------------------------------
# region notes
# Zone C = Paris / Île-de-France. Affects day-type dispatch (JOHV/JOVS/WKD).
# Source: Ministère de l'Éducation Nationale, calendrier scolaire 2023-2024.
# Extend _ZONE_C_HOLIDAYS and _PUBLIC_HOLIDAYS for other simulation years.
# endregion

_ZONE_C_HOLIDAYS_2024 = [
    ("2024-02-17", "2024-03-04"),
    ("2024-04-13", "2024-04-29"),
    ("2024-07-06", "2024-09-02"),
    ("2024-10-19", "2024-11-04"),
    ("2024-12-21", "2025-01-06"),
]

_PUBLIC_HOLIDAYS_2024 = [
    "2024-01-01", "2024-04-01", "2024-05-01", "2024-05-08",
    "2024-05-09", "2024-05-20", "2024-07-14", "2024-08-15",
    "2024-11-01", "2024-11-11", "2024-12-25",
]


def _is_holiday(ts: pd.Timestamp) -> bool:
    date_str = ts.strftime("%Y-%m-%d")
    if date_str in _PUBLIC_HOLIDAYS_2024:
        return True
    for start, end in _ZONE_C_HOLIDAYS_2024:
        if start <= date_str <= end:
            return True
    return False


def _day_type(ts: pd.Timestamp) -> str:
    """JOHV / JOVS / WKD dispatch."""
    if ts.dayofweek >= 5 or _is_holiday(ts):
        return "WKD"
    date_str = ts.strftime("%Y-%m-%d")
    for start, end in _ZONE_C_HOLIDAYS_2024:
        if start <= date_str <= end:
            return "JOVS"
    return "JOHV"


# -----------------------------------------------------------------------------
# 2. LOAD OCCUPANCY PROFILES
# -----------------------------------------------------------------------------

def load_profiles(path: str) -> dict:
    """Read Defence_Occupation_Normalised.xlsx → dict[str, np.ndarray(24)], values 0–1."""
    df = pd.read_excel(path)
    return {key: df[key].values[:24].astype(float) for key in ("JOHV", "JOVS", "WKD")}


# -----------------------------------------------------------------------------
# 3. BUILD Q AND N ARRAYS
# -----------------------------------------------------------------------------

def build_Q_array(
    dates: pd.DatetimeIndex,
    profiles: dict,
    watts_per_person: float = WATTS_SENSIBLE_PP,
    baseline_w: float = BASELINE_W,
    people_peak: int = PEOPLE_PEAK,
) -> tuple:
    """
    Returns (Q_array [W], n_people_array [persons]) for a datetime index.
    Q = baseline + n_people * watts_per_person.
    """
    n = len(dates)
    Q_array        = np.empty(n)
    n_people_array = np.empty(n)
    for i, ts in enumerate(dates):
        frac             = profiles[_day_type(ts)][ts.hour]
        n_people         = frac * people_peak
        Q_array[i]       = baseline_w + n_people * watts_per_person
        n_people_array[i] = n_people
    return Q_array, n_people_array


# -----------------------------------------------------------------------------
# 4. INFILTRATION RATE  V_inf [m³/s]
# -----------------------------------------------------------------------------
# region notes
# V_inf = V_cycle / T_headway = (eta * V_platform) / T_headway
# Hour boundaries are engineering estimates — DUP does not publish hourly headways.
# WKD: off-peak all day. JOHV/JOVS: peak at 06-09h and 17-20h.
# Night (01h–05h): no trains, V_inf = 0. Hours 23–00h kept as offpeak (trains thinning
# out — real headway ~15 min, modeled as 4 min offpeak; noted as limitation in README).
# endregion

_WEEKDAY_REGIME = {
    0:  "offpeak",
    **{h: "night"   for h in range(1, 5)},
    5:  "offpeak",
    **{h: "peak"    for h in range(6, 9)},
    **{h: "offpeak" for h in range(9, 17)},
    **{h: "peak"    for h in range(17, 20)},
    **{h: "offpeak" for h in range(20, 24)},
}

_WKD_REGIME = {h: ("night" if 1 <= h < 5 else "offpeak") for h in range(24)}


def v_inf_m3s(hour: int, day_type: str) -> float:
    """Time-averaged infiltration rate [m³/s] for given hour and day_type."""
    regime = _WKD_REGIME[hour] if day_type == "WKD" else _WEEKDAY_REGIME[hour]
    if regime == "night":
        return 0.0
    elif regime == "peak":
        return V_CYCLE_M3 / HEADWAY_PEAK_S       # 138.6 / 120 = 1.155 m³/s
    else:
        return V_CYCLE_M3 / HEADWAY_OFFPEAK_S    # 138.6 / 240 = 0.578 m³/s


# -----------------------------------------------------------------------------
# 5. SMOKE TEST
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    import os
    print("V_inf checks:")
    print(f"  Peak 08h JOHV:  {v_inf_m3s(8,  'JOHV'):.3f} m³/s")
    print(f"  Offpk 11h JOHV: {v_inf_m3s(11, 'JOHV'):.3f} m³/s")
    print(f"  Night 03h JOHV: {v_inf_m3s(3,  'JOHV'):.3f} m³/s")
    print(f"  WKD   08h:      {v_inf_m3s(8,  'WKD'):.3f} m³/s")

    xlsx_path = os.path.join("data", "raw", "Defense_Occupation_Normalised.xlsx")
    if os.path.exists(xlsx_path):
        profiles = load_profiles(xlsx_path)
        dates = pd.date_range("2024-07-01", periods=168, freq="h")
        Q, n = build_Q_array(dates, profiles)
        print(f"Q: {Q.min()/1000:.1f}–{Q.max()/1000:.1f} kW | n: {n.min():.0f}–{n.max():.0f} p")