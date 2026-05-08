# =============================================================================
# occupancy.py — Occupancy profiles and infiltration rate for one platform zone
#
# Provides:
#   load_profiles(path)            → dict of hourly occupancy profiles (0–1)
#   build_Q_array(dates, profiles) → (Q_array, n_people_array) for a datetime index
#   v_inf_m3s(hour, day_type)      → time-averaged infiltration rate [m³/s]
#
# Scope: one platform zone (one side, one AHU). 250-person peak.
# Day-type dispatch: JOHV / JOVS / WKD with Zone C 2024 holiday calendar.
# =============================================================================

import numpy as np
import pandas as pd

from constants import (
    PEOPLE_PEAK,
    WATTS_SENSIBLE_PP,
    BASELINE_W,
    ETA_INF,
    V_CYCLE_M3,
    HEADWAY_PEAK_S,
    HEADWAY_OFFPEAK_S,
    HEADWAY_NIGHT_S,
)


# -----------------------------------------------------------------------------
# SECTION 1 — Zone C 2024 holiday calendar (hardcoded)
# -----------------------------------------------------------------------------
# French school holiday zones affect ridership patterns on the RATP/GPE network.
# Zone C covers Paris and Île-de-France. Holidays shift weekday profiles toward
# a reduced pattern (similar to weekend behaviour for commuter-dominated stations).
# Source: Ministère de l'Éducation Nationale, calendrier scolaire 2023-2024:
# https://www.education.gouv.fr/les-dates-des-vacances-scolaires-100660
# Only 2024 dates are hardcoded — extend for other simulation years.

_ZONE_C_HOLIDAYS_2024 = [
    # Vacances d'hiver
    ("2024-02-17", "2024-03-04"),
    # Vacances de printemps
    ("2024-04-13", "2024-04-29"),
    # Vacances d'été
    ("2024-07-06", "2024-09-02"),
    # Vacances de la Toussaint
    ("2024-10-19", "2024-11-04"),
    # Vacances de Noël
    ("2024-12-21", "2025-01-06"),
]

_PUBLIC_HOLIDAYS_2024 = [
    "2024-01-01",  # Jour de l'An
    "2024-04-01",  # Lundi de Pâques
    "2024-05-01",  # Fête du Travail
    "2024-05-08",  # Victoire 1945
    "2024-05-09",  # Ascension
    "2024-05-20",  # Lundi de Pentecôte
    "2024-07-14",  # Fête Nationale
    "2024-08-15",  # Assomption
    "2024-11-01",  # Toussaint
    "2024-11-11",  # Armistice
    "2024-12-25",  # Noël
]


def _is_holiday(ts: pd.Timestamp) -> bool:
    """Return True if ts falls in a Zone C school holiday or public holiday."""
    date_str = ts.strftime("%Y-%m-%d")
    if date_str in _PUBLIC_HOLIDAYS_2024:
        return True
    for start, end in _ZONE_C_HOLIDAYS_2024:
        if start <= date_str <= end:
            return True
    return False


def _day_type(ts: pd.Timestamp) -> str:
    """
    Map a timestamp to its RATP day-type code:
      JOHV — Jour Ouvrable Hors Vacances (standard weekday, school in session)
      JOVS — Jour Ouvrable Vacances Scolaires (weekday, school holiday)
      WKD  — Weekend / public holiday
    """
    if ts.dayofweek >= 5 or _is_holiday(ts):
        return "WKD"
    # Weekday during school holiday period → JOVS
    date_str = ts.strftime("%Y-%m-%d")
    for start, end in _ZONE_C_HOLIDAYS_2024:
        if start <= date_str <= end:
            return "JOVS"
    return "JOHV"


# -----------------------------------------------------------------------------
# SECTION 2 — Load normalised occupancy profiles from Excel
# -----------------------------------------------------------------------------

def load_profiles(path: str) -> dict:
    """
    Read the normalised RATP hourly occupancy profiles.

    Parameters
    ----------
    path : str
        Path to Defense_Occupation_Normalised.xlsx.

    Returns
    -------
    dict[str, np.ndarray(24)]
        Keys: 'JOHV', 'JOVS', 'WKD'.
        Values: hourly fraction 0–1, where 1.0 = JOHV 18h peak (250 persons).
    """
    df = pd.read_excel(path)  # columns: heure, JOHV, JOVS, WKD
    profiles = {}
    for key in ("JOHV", "JOVS", "WKD"):
        profiles[key] = df[key].values[:24].astype(float)
        # already normalized against JOHV 18h peak in the Excel — do not rescale
    return profiles


# -----------------------------------------------------------------------------
# SECTION 3 — Build Q_array and n_people_array from a datetime index
# -----------------------------------------------------------------------------

def build_Q_array(
    dates: pd.DatetimeIndex,
    profiles: dict,
    watts_per_person: float = WATTS_SENSIBLE_PP,
    baseline_w: float = BASELINE_W,
    people_peak: int = PEOPLE_PEAK,
) -> tuple:
    """
    Build hourly internal heat gain and headcount arrays.

    Parameters
    ----------
    dates : pd.DatetimeIndex
        Hourly timestamps for the simulation window.
    profiles : dict
        Output of load_profiles().
    watts_per_person : float
        Sensible heat per occupant [W]. Default from constants.
    baseline_w : float
        Continuous equipment load [W]. Default from constants.
    people_peak : int
        Peak headcount for the zone. Default from constants.

    Returns
    -------
    Q_array : np.ndarray [W]
        Total internal heat gain at each timestep.
    n_people_array : np.ndarray [persons]
        Headcount at each timestep.
    """
    n = len(dates)
    Q_array       = np.empty(n)
    n_people_array = np.empty(n)

    for i, ts in enumerate(dates):
        dtype   = _day_type(ts)
        frac    = profiles[dtype][ts.hour]
        n_people = frac * people_peak
        Q_array[i]        = baseline_w + n_people * watts_per_person
        n_people_array[i] = n_people

    return Q_array, n_people_array


# -----------------------------------------------------------------------------
# SECTION 4 — Time-averaged infiltration rate V̇_inf [m³/s]
# -----------------------------------------------------------------------------
# Physical basis: see constants.py §6 for full derivation and justification.
#
# V̇_inf = V_cycle / T_headway
#        = (η × V_platform) / T_headway
#
# Train frequency varies by hour and day-type. Three regimes are defined:
#
#   PEAK hours     → HEADWAY_PEAK_S    (2 min, DUP 2012)
#   OFF-PEAK hours → HEADWAY_OFFPEAK_S (4 min, [ASSUMPTION])
#   NIGHT hours    → no service → V̇_inf = 0 + residual gap infiltration
#
# Dispatch by hour (applies to all day-types — frequency profile is the
# same shape; absolute headway differs by day-type only at peak):
#
#   Night  (00–05):  no service
#   Ramp   (05–06):  off-peak headway (first trains)
#   Peak   (06–09):  peak headway
#   Normal (09–17):  off-peak headway
#   Peak   (17–20):  peak headway
#   Normal (20–22):  off-peak headway
#   Night  (22–24):  last trains / no service → off-peak, then zero
#
# [ASSUMPTION] — these hour boundaries are engineering estimates consistent
# with Île-de-France Mobilités published service patterns for GPE.
# Exact hourly headways are not publicly detailed in the DUP.
# WKD headways are assumed to follow the off-peak profile throughout the day
# (no WKD peak headway published for ligne 15 in public documents).
#
# Residual gap infiltration at night (PSDs closed, no trains):
# [ASSUMPTION] — treated as zero in this model for simplicity.
# Real residual infiltration through PSD gaps exists but is small and
# not quantifiable without PSD leakage area data (a manufacturer spec
# not publicly available). Flagged for future refinement.
# [SOBOL — η range 0.10–0.25 covers the dominant uncertainty]

# Hour → regime lookup for weekdays
_HOUR_TO_REGIME_WEEKDAY = {
    0: "night",  1: "night",  2: "night",  3: "night",  4: "night",
    5: "offpeak",
    6: "peak",   7: "peak",   8: "peak",   9: "offpeak",
    10: "offpeak", 11: "offpeak", 12: "offpeak", 13: "offpeak",
    14: "offpeak", 15: "offpeak", 16: "offpeak",
    17: "peak",  18: "peak",  19: "peak",
    20: "offpeak", 21: "offpeak",
    22: "offpeak", 23: "night",
}

# Weekend / holiday: off-peak all day, no peak service
_HOUR_TO_REGIME_WKD = {h: ("offpeak" if 6 <= h <= 22 else "night")
                        for h in range(24)}


def v_inf_m3s(hour: int, day_type: str) -> float:
    """
    Time-averaged platform infiltration rate at a given hour and day-type.

    Parameters
    ----------
    hour : int
        Hour of day (0–23).
    day_type : str
        'JOHV', 'JOVS', or 'WKD'.

    Returns
    -------
    float
        V̇_inf in m³/s (time-averaged over the full hour).

    Notes
    -----
    JOHV and JOVS follow the weekday regime table (peak hours during
    AM and PM rush). WKD uses a flat off-peak profile throughout.
    The return value feeds directly into the ODE slope function as:
        Q_inf = RHO_CP_AIR * v_inf * (T_tun - T_in)
    """
    if day_type == "WKD":
        regime = _HOUR_TO_REGIME_WKD[hour]
    else:
        # JOHV and JOVS share the same hour pattern; frequency is the same
        # (DUP does not distinguish headway between JOHV and JOVS).
        regime = _HOUR_TO_REGIME_WEEKDAY[hour]

    if regime == "night":
        return 0.0
    elif regime == "peak":
        return V_CYCLE_M3 / HEADWAY_PEAK_S        # 138.6 / 120 = 1.155 m³/s
    else:  # offpeak
        return V_CYCLE_M3 / HEADWAY_OFFPEAK_S     # 138.6 / 240 = 0.578 m³/s


# -----------------------------------------------------------------------------
# SECTION 5 — Smoke test
# -----------------------------------------------------------------------------
if __name__ == "__main__":
    import os

    print("=== occupancy.py smoke test ===\n")

    # --- v_inf_m3s ---
    print("V̇_inf by regime:")
    print(f"  Weekday peak  (08h JOHV): {v_inf_m3s(8,  'JOHV'):.3f} m³/s")
    print(f"  Weekday offpk (11h JOHV): {v_inf_m3s(11, 'JOHV'):.3f} m³/s")
    print(f"  Night         (03h JOHV): {v_inf_m3s(3,  'JOHV'):.3f} m³/s")
    print(f"  Weekend peak  (08h WKD):  {v_inf_m3s(8,  'WKD'):.3f} m³/s")
    print()

    # --- load_profiles ---
    xlsx_path = os.path.join("data", "raw", "Defense_Occupation_Normalised.xlsx")
    if os.path.exists(xlsx_path):
        profiles = load_profiles(xlsx_path)
        print("Profiles loaded:")
        for k, v in profiles.items():
            print(f"  {k}: min={v.min():.2f}, max={v.max():.2f}, peak_hour={v.argmax()}h")

        # --- build_Q_array ---
        dates = pd.date_range("2024-07-01", periods=168, freq="h")
        Q, n = build_Q_array(dates, profiles)
        print(f"\nbuild_Q_array (July 2024, 1 week):")
        print(f"  Q range: {Q.min()/1000:.1f}–{Q.max()/1000:.1f} kW")
        print(f"  n range: {n.min():.0f}–{n.max():.0f} persons")
    else:
        print(f"[SKIP] Excel not found at {xlsx_path} — skipping profile test.")
