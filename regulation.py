# =============================================================================
# regulation.py — HVAC regulation logic for one platform zone
#
# Contains:
#   T_setpoint(T_ext)                     → target indoor temperature [°C]
#   airflow_total(n_people)               → total AHU airflow [m³/h]
#   dT_dt(t, T, ...)                      → ODE slope function [K/s]
#   build_Q_hvac_array(T_in, T_ext, ...)  → post-hoc HVAC decomposition
#
# ODE (updated Session 11):
#   C · dT_in/dt = (UA_f + ρcp·V̇_inf)·(T_tun − T_in)
#                + UA_s·(T_soil − T_in)
#                + Q_int
#                − Q_hvac
#
# T_ext no longer appears as a direct envelope boundary.
# T_tun = T_ext + T_TUN_OFFSET_C  (computed inside dT_dt)
# T_soil = T_SOIL_C               (constant, from constants.py)
# V̇_inf = v_inf_m3s(hour, day_type) (from occupancy.py, dynamic)
#
# Scope: one platform zone (one side, one AHU). 250-person peak.
# =============================================================================

import numpy as np
import pandas as pd

from constants import (
    # Envelope
    UA_FACADE_W_K,
    UA_SOIL_W_K,
    RHO_CP_AIR_J_M3_K,
    # Boundary temperatures
    T_SOIL_C,
    T_TUN_OFFSET_C,
    # Capacitance
    C_TOTAL_J_K,
    # Airflow
    AIRFLOW_MIN_M3H,
    AIRFLOW_MAX_M3H,
    AIRFLOW_PER_PERSON_M3H,
    AIRFLOW_OVERPRESSURE_M3H,
    # Regulation setpoint boundaries
    T_ANTIFREEZE_C,
    T_HEAT_LINEAR_LOW_C,
    T_HEAT_LINEAR_HIGH_C,
    T_HEAT_FIXED_C,
    T_DEAD_LOW_C,
    T_DEAD_HIGH_C,
    T_COOL_FIXED_C,
    T_COOL_LINEAR_HIGH_C,
    T_BLOW_HEAT_C,
    T_BLOW_COOL_C,
    # Occupancy
    PEOPLE_PEAK,
)
from occupancy import v_inf_m3s, _day_type


# -----------------------------------------------------------------------------
# SECTION 1 — Setpoint law: T_set = f(T_ext)
# -----------------------------------------------------------------------------
# Five-zone piecewise law driven by outdoor temperature T_ext.
# T_ext is still the control input for the setpoint law — the controller
# reads outdoor air temperature to decide what indoor temperature to target.
# (T_ext is NOT a thermal boundary condition on the envelope anymore, but it
# remains the scheduling variable for the regulation strategy.)
#
# Zone definitions:
#   T_ext < -1°C     → Anti-freeze: T_set = 5°C (hard floor)
#   -1 ≤ T_ext < 6   → Linear heating offset: T_set = T_ext + 6
#   6 ≤ T_ext < 12   → Fixed heating target: T_set = 12°C
#   12 ≤ T_ext ≤ 26  → Dead band: T_set = NaN (no active thermal target)
#   26 < T_ext ≤ 31  → Fixed cooling target: T_set = 26°C
#   T_ext > 31       → Linear cooling offset: T_set = T_ext − 5

def T_setpoint(T_ext: float) -> float:
    """
    Compute the indoor air temperature setpoint from outdoor temperature.

    Parameters
    ----------
    T_ext : float
        Outdoor air temperature [°C].

    Returns
    -------
    float or np.nan
        Target indoor temperature [°C].
        Returns np.nan in the dead band (no active thermal control).
    """
    if T_ext < T_HEAT_LINEAR_LOW_C:
        return T_ANTIFREEZE_C                          # anti-freeze floor
    elif T_ext < T_HEAT_LINEAR_HIGH_C:
        return T_ext + 6.0                             # linear heating offset
    elif T_ext < T_DEAD_LOW_C:
        return T_HEAT_FIXED_C                          # fixed heating target
    elif T_ext <= T_DEAD_HIGH_C:
        return np.nan                                  # dead band
    elif T_ext <= T_COOL_LINEAR_HIGH_C:
        return T_COOL_FIXED_C                          # fixed cooling target
    else:
        return T_ext - 5.0                             # linear cooling offset


# -----------------------------------------------------------------------------
# SECTION 2 — Airflow modulation: Q_air = f(n_people)
# -----------------------------------------------------------------------------
# Airflow is proportional to instantaneous occupancy, bounded between:
#   AIRFLOW_MIN: overpressure floor only (zero occupancy)
#   AIRFLOW_MAX: peak demand + safety margin
#
# Additive sizing: Q_total = Q_overpressure + Q_occupancy
# See constants.py §7-8 and parameters.md §7-8 for full justification.

def airflow_total(n_people: float) -> float:
    """
    Compute total AHU airflow for a given instantaneous headcount.

    Parameters
    ----------
    n_people : float
        Current number of persons in the zone.

    Returns
    -------
    float
        Total airflow [m³/h], bounded [AIRFLOW_MIN, AIRFLOW_MAX].
    """
    n_clamped = min(n_people, PEOPLE_PEAK)
    Q = AIRFLOW_OVERPRESSURE_M3H + n_clamped * AIRFLOW_PER_PERSON_M3H
    return float(np.clip(Q, AIRFLOW_MIN_M3H, AIRFLOW_MAX_M3H))


# -----------------------------------------------------------------------------
# SECTION 3 — ODE slope function: dT_in/dt
# -----------------------------------------------------------------------------
# Full equation (Session 11 update):
#
#   C · dT/dt = (UA_f + ρcp·V̇_inf)·(T_tun − T)
#             +  UA_s·(T_soil − T)
#             +  Q_int
#             −  Q_hvac
#
# Where:
#   T_tun  = T_ext + T_TUN_OFFSET_C        (tunnel air temperature)
#   T_soil = T_SOIL_C                       (constant ground temperature)
#   V̇_inf  = v_inf_m3s(hour, day_type)     (dynamic, from occupancy.py)
#   Q_hvac = ρcp · Q_air_m3s · ΔT_blow     (AHU heating/cooling power)
#   ΔT_blow = clip(T_set − T_in, −12, +5)  (asymmetric saturation)
#
# ΔT_blow saturation:
#   +5 K max (cooling): T_ext > 31, T_set = T_ext − 5 → ΔT = 5 by construction
#   −12 K max (heating): T_ext = −7, T_set = 5 → ΔT = 5 − (−7) = 12
#   Clip ensures physical plausibility without the water regime cap (pending S11).
#
# Dead band behaviour:
#   T_set = NaN → Q_hvac = 0. But ventilation continues at minimum airflow.
#   The AHU blows outdoor air (conditioned to T_ext approximately) into the zone.
#   This ventilation term is captured in the envelope equation:
#   the AHU draws outdoor air at T_ext, which reaches the zone at ~T_ext after
#   the duct. Net effect: Q_vent = ρcp · Q_air_m3s · (T_ext − T_in).
#   This term is included explicitly in dT_dt via the Q_hvac branch for dead band.

def dT_dt(
    t: float,
    T: list,
    t_array: np.ndarray,
    T_ext_array: np.ndarray,
    Q_int_array: np.ndarray,
    n_people_array: np.ndarray,
    dates: pd.DatetimeIndex,
) -> list:
    """
    ODE slope function for the lumped-capacitance platform thermal model.

    Parameters
    ----------
    t : float
        Current solver time [s], measured from t_array[0].
    T : list[float]
        Current state: [T_in] in °C.
    t_array : np.ndarray
        Time array [s] aligned with the simulation window.
    T_ext_array : np.ndarray
        Outdoor air temperature [°C] at each timestep.
    Q_int_array : np.ndarray
        Internal heat gain [W] at each timestep (occupancy + equipment).
    n_people_array : np.ndarray
        Headcount [persons] at each timestep.
    dates : pd.DatetimeIndex
        Datetime index corresponding to t_array (for day-type dispatch).

    Returns
    -------
    list[float]
        [dT_in/dt] in K/s.
    """
    T_in = float(T[0])

    # --- Interpolate inputs at current solver time ---
    T_ext   = float(np.interp(t, t_array, T_ext_array))
    Q_int   = float(np.interp(t, t_array, Q_int_array))
    n_ppl   = float(np.interp(t, t_array, n_people_array))

    # --- Identify current timestep index for day-type dispatch ---
    idx     = int(np.clip(np.searchsorted(t_array, t), 0, len(dates) - 1))
    ts      = dates[idx]
    dtype   = _day_type(ts)
    hour    = ts.hour

    # --- Boundary temperatures ---
    T_tun  = T_ext + T_TUN_OFFSET_C     # tunnel air [°C]
    T_soil = T_SOIL_C                   # stable ground temperature [°C]

    # --- Infiltration rate ---
    V_inf  = v_inf_m3s(hour, dtype)     # m³/s, time-averaged

    # --- Effective facade conductance (conduction + infiltration) ---
    UA_eff_tun = UA_FACADE_W_K + RHO_CP_AIR_J_M3_K * V_inf    # W/K

    # --- Envelope heat transfer ---
    Q_facade = UA_eff_tun  * (T_tun  - T_in)   # W (positive = heat gain)
    Q_soil   = UA_SOIL_W_K * (T_soil - T_in)   # W (positive = heat gain)

    # --- HVAC power ---
    T_set = T_setpoint(T_ext)
    Q_air_m3h = airflow_total(n_ppl)
    Q_air_m3s = Q_air_m3h / 3600.0

    if np.isnan(T_set):
        # Dead band: no active thermal control.
        # AHU blows outdoor air at ~T_ext into the zone.
        # Net ventilation heat exchange: ρcp · Q_air · (T_ext − T_in)
        Q_hvac = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_ext)
    else:
    # Fixed T_blow: Q_hvac = rho_cp * Q_air * (T_in - T_blow)
    # Positive = cooling (T_in > T_blow), negative = heating (T_in < T_blow)
    # Controller shuts off if T_in already on the right side of T_set
        if T_set > T_in:   # heating needed
            T_blow = T_BLOW_HEAT_C
            Q_hvac = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_blow)
            if T_in >= T_set:  # already warm enough, shut off
                Q_hvac = 0.0
        else:              # cooling needed
            T_blow = T_BLOW_COOL_C
            Q_hvac = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_blow)
            if T_in <= T_set:  # already cool enough, shut off
                Q_hvac = 0.0

    # --- ODE ---
    dTdt = (Q_facade + Q_soil + Q_int - Q_hvac) / C_TOTAL_J_K

    return [dTdt]


# -----------------------------------------------------------------------------
# SECTION 4 — Post-hoc HVAC decomposition (for plotting only)
# -----------------------------------------------------------------------------
# Reconstructs Q_heat, Q_cool, Q_vent from stored T_in and T_ext arrays.
# DO NOT call from inside dT_dt — that would be O(N²).
# Call once after solve_ivp returns, on the full output arrays.

def build_Q_hvac_array(
    T_in_array: np.ndarray,
    T_ext_array: np.ndarray,
    n_people_array: np.ndarray,
    dates: pd.DatetimeIndex,
) -> tuple:
    """
    Decompose HVAC power into heating, cooling, and ventilation components.

    Parameters
    ----------
    T_in_array : np.ndarray
    T_ext_array : np.ndarray
    n_people_array : np.ndarray
    dates : pd.DatetimeIndex

    Returns
    -------
    Q_hvac_total, Q_heat, Q_cool, Q_vent : np.ndarray [W]
        Mutually exclusive. Q_heat < 0 (adds heat), Q_cool > 0 (removes heat).
    """
    n = len(T_in_array)
    Q_total = np.zeros(n)
    Q_heat  = np.zeros(n)
    Q_cool  = np.zeros(n)
    Q_vent  = np.zeros(n)

    for i, (T_in, T_ext, n_ppl, ts) in enumerate(
            zip(T_in_array, T_ext_array, n_people_array, dates)):

        T_set     = T_setpoint(T_ext)
        Q_air_m3h = airflow_total(n_ppl)
        Q_air_m3s = Q_air_m3h / 3600.0

        if np.isnan(T_set):
            # Dead band: ventilation only, supply air ≈ T_ext
            q = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_ext)
            Q_vent[i]  = q
            Q_total[i] = q
        else:
            if T_set > T_in:
                # Heating needed: blow hot air
                T_blow = T_BLOW_HEAT_C
                q = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_blow)  # negative
                q = q if T_in < T_set else 0.0   # shut off if already warm enough
                Q_heat[i]  = q
            else:
                # Cooling needed: blow cold air
                T_blow = T_BLOW_COOL_C
                q = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_blow)  # positive
                q = q if T_in > T_set else 0.0   # shut off if already cool enough
                Q_cool[i]  = q
            Q_total[i] = q

    return Q_total, Q_heat, Q_cool, Q_vent