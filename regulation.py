# =============================================================================
# regulation.py — HVAC regulation logic, one platform zone
# Functions: T_setpoint, airflow_total, dT_dt, T_hot_water_supply,
#            T_cold_water_supply, build_Q_hvac_array
# =============================================================================

import numpy as np
import pandas as pd

from constants import (
    UA_FACADE_W_K, UA_TUN_WALL_W_K, UA_SOIL_W_K, RHO_CP_AIR_J_M3_K,
    T_SOIL_C, T_TUN_OFFSET_DAY_C, T_TUN_OFFSET_NIGHT_C,
    C_TOTAL_J_K,
    AIRFLOW_MIN_M3H, AIRFLOW_MAX_M3H, AIRFLOW_PER_PERSON_M3H, AIRFLOW_OVERPRESSURE_M3H,
    Q_STAIR_M3S, T_STAIR_COLD_C, F_STAIR_COLD, F_STAIR_NIGHT,
    Q_CURTAIN_AIR_M3S, DT_JET_K, ETA_IN_CURTAIN,
    T_EXT_HEAT_LOW_C, T_HEAT_LOW_C, T_HEAT_HIGH_C,
    T_DEAD_LOW_C, T_DEAD_HIGH_C,
    T_COOL_FIXED_C, T_COOL_DELTA_C, T_EXT_DELTA_C,
    T_BLOW_HEAT_C, T_BLOW_COOL_C, T_NIGHT_SETBACK_C,
    PEOPLE_PEAK,
    T_HW_EXT_LOW_C, T_HW_EXT_HIGH_C, T_HW_EXT_HYST_C, T_HW_SUPPLY_MAX, T_HW_SUPPLY_MIN,
    T_CW_EXT_LOW_C, T_CW_EXT_HYST_C, T_CW_SUPPLY_C,
    FRAC_RETURN_AIR, CP_GLYCOL_J_KG_K, RHO_GLYCOL_KG_M3, DT_WATER_HEAT_K, DT_WATER_COOL_K,
)
from occupancy import v_inf_m3s, _day_type


# -----------------------------------------------------------------------------
# 1. SETPOINT LAW  T_set = f(T_ext, hour)
# -----------------------------------------------------------------------------
# region notes
# Night (01h–05h): anti-freeze only → T_set = 5°C if heating needed, else NaN.
# Service hours, 5-zone outdoor compensation:
#   T_ext ≤  5°C       → 18°C
#   5 < T_ext < 15°C   → 18–20°C (linear)
#   15 ≤ T_ext ≤ 22°C  → dead band (NaN)
#   22 < T_ext ≤ 32°C  → 26°C
#   T_ext > 32°C       → min(27, T_ext − 6)
# endregion

def T_setpoint(T_ext: float, hour: int = 12) -> float:
    """T_set [°C] from T_ext [°C] and hour. Returns np.nan in dead band.
    Night (01h–05h): anti-freeze at 5°C only."""
    is_night = 1 <= hour < 5

    if T_ext <= T_EXT_HEAT_LOW_C:
        return T_NIGHT_SETBACK_C if is_night else T_HEAT_LOW_C
    elif T_ext < T_DEAD_LOW_C:
        if is_night:
            return T_NIGHT_SETBACK_C
        return float(np.interp(T_ext,
            [T_EXT_HEAT_LOW_C, T_DEAD_LOW_C],
            [T_HEAT_LOW_C, T_HEAT_HIGH_C]))        # 18→20°C linear
    elif T_ext <= T_DEAD_HIGH_C:
        return np.nan                               # dead band
    elif T_ext <= T_EXT_DELTA_C:
        return T_COOL_FIXED_C                       # 26°C
    else:
        return T_ext - T_COOL_DELTA_C             # T_ext − 6°C

# -----------------------------------------------------------------------------
# 1b. STAIRCASE FLOW MODULATION  q_stair = f(T_ext, hour)
# -----------------------------------------------------------------------------

def q_stair_m3s(T_ext: float, hour: int) -> float:
    """Effective stair airflow [m³/s] after curtain / shutter modulation."""
    if 1 <= hour < 5:
        return F_STAIR_NIGHT * Q_STAIR_M3S          # metal door closed
    elif T_ext < T_STAIR_COLD_C:
        return F_STAIR_COLD * Q_STAIR_M3S            # air curtain active
    else:
        return Q_STAIR_M3S                            # fully open


def curtain_active(T_ext: float, hour: int) -> bool:
    """True when air curtain is running (service hours + cold weather)."""
    return not (1 <= hour < 5) and T_ext < T_STAIR_COLD_C


# -----------------------------------------------------------------------------
# 2. AIRFLOW  Q_air = f(n_people)
# -----------------------------------------------------------------------------

def airflow_total(n_people: float) -> float:
    """Total AHU airflow [m³/h] for current headcount."""
    n_clamped = min(n_people, PEOPLE_PEAK)
    Q = AIRFLOW_OVERPRESSURE_M3H + n_clamped * AIRFLOW_PER_PERSON_M3H
    return float(np.clip(Q, AIRFLOW_MIN_M3H, AIRFLOW_MAX_M3H))


# -----------------------------------------------------------------------------
# 3. ODE  dT_in/dt
# -----------------------------------------------------------------------------
# region equation
# C·dT/dt = (UA_facade + UA_tun_wall + ρcp·V_inf)·(T_tun - T)
#          + UA_soil·(T_soil - T)
#          + Q_int
#          + ρcp·Q_stair(T_ext,hour)·(T_ext - T)
#          + Q_curtain_zone   (when curtain active)
#          - Q_hvac
#
# T_tun = T_ext + offset (10°C service / 5°C night 01h–05h)
# Q_hvac positive = heat removed from zone
# Dead band: pure ventilation (AHU blows outdoor air)
# Heating/cooling airflow boost: if hygiene flow can't meet load, boost up to AIRFLOW_MAX
# endregion

def dT_dt(
    t: float,
    T: list,
    t_array: np.ndarray,
    T_ext_array: np.ndarray,
    Q_int_array: np.ndarray,
    n_people_array: np.ndarray,
    dates: pd.DatetimeIndex,
    water_state: dict,
) -> list:
    """ODE slope [K/s]. water_state = {'heating': bool, 'cooling': bool} — mutable, persists."""
    T_in = float(T[0])

    T_ext = float(np.interp(t, t_array, T_ext_array))
    Q_int = float(np.interp(t, t_array, Q_int_array))
    n_ppl = float(np.interp(t, t_array, n_people_array))

    idx   = int(np.clip(np.searchsorted(t_array, t), 0, len(dates) - 1))
    ts    = dates[idx]
    dtype = _day_type(ts)
    hour  = ts.hour

    # Boundary temperatures — night = 01h–05h
    is_night = 1 <= hour < 5
    offset   = T_TUN_OFFSET_NIGHT_C if is_night else T_TUN_OFFSET_DAY_C
    T_tun    = T_ext + offset
    T_soil   = T_SOIL_C

    # Envelope gains (W)
    V_inf      = v_inf_m3s(hour, dtype)
    UA_eff_tun = UA_FACADE_W_K + UA_TUN_WALL_W_K + RHO_CP_AIR_J_M3_K * V_inf
    Q_facade   = UA_eff_tun  * (T_tun  - T_in)
    Q_soil     = UA_SOIL_W_K * (T_soil - T_in)

    # Staircase — modulated by curtain/shutter
    Q_stair_flow = q_stair_m3s(T_ext, hour)
    Q_stair      = RHO_CP_AIR_J_M3_K * Q_stair_flow * (T_ext - T_in)

    # Air curtain heat spill into zone
    Q_curtain_zone = 0.0
    if curtain_active(T_ext, hour):
        P_heat_curtain = RHO_CP_AIR_J_M3_K * Q_CURTAIN_AIR_M3S * DT_JET_K
        Q_curtain_zone = ETA_IN_CURTAIN * P_heat_curtain

    # HVAC
    T_set     = T_setpoint(T_ext, hour)
    Q_air_m3s = airflow_total(n_ppl) / 3600.0

    T_hw = T_hot_water_supply(T_ext, water_state["heating"])
    T_cw = T_cold_water_supply(T_ext, water_state["cooling"])
    water_state["heating"] = T_hw is not None
    water_state["cooling"] = T_cw is not None

    if np.isnan(T_set):
        # Dead band — pure ventilation with outdoor air
        Q_hvac = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_ext) * (1 - FRAC_RETURN_AIR)

    elif T_set > T_in:
        # Heating — airflow boost if hygiene flow can't keep up
        if water_state["heating"] and T_in < T_set:
            Q_gains = Q_facade + Q_soil + Q_int + Q_stair + Q_curtain_zone  # negative in winter
            heating_capacity = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_BLOW_HEAT_C - T_in)  # positive
            if Q_gains + heating_capacity < 0:
                # Hygiene flow can't cover losses — boost
                Q_air_required = -Q_gains / (RHO_CP_AIR_J_M3_K * (T_BLOW_HEAT_C - T_in))
                Q_air_m3s = min(Q_air_required, AIRFLOW_MAX_M3H / 3600.0)
            Q_hvac = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_BLOW_HEAT_C) # return air already accounted for in T_blow_heat_C
        else:
            Q_hvac = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_ext) * (1 - FRAC_RETURN_AIR)

    else:
        # Cooling — airflow boost if hygiene flow can't meet load
        if water_state["cooling"] and T_in > T_BLOW_COOL_C:
            Q_load = Q_facade + Q_soil + Q_int + Q_stair + Q_curtain_zone
            Q_cool_hygiene = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_BLOW_COOL_C)
            if Q_load > Q_cool_hygiene:
                Q_air_required = Q_load / (RHO_CP_AIR_J_M3_K * (T_in - T_BLOW_COOL_C))
                Q_air_m3s = min(Q_air_required, AIRFLOW_MAX_M3H / 3600.0)
            Q_hvac = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_BLOW_COOL_C) # return air already accounted for in T_blow_cool_C
        else:
            Q_hvac = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_ext) * (1 - FRAC_RETURN_AIR)

    dTdt = (Q_facade + Q_soil + Q_int + Q_stair + Q_curtain_zone - Q_hvac) / C_TOTAL_J_K
    return [dTdt]


# -----------------------------------------------------------------------------
# 4. WATER REGIME  T_water = f(T_ext)
# -----------------------------------------------------------------------------
# region notes
# Hot:  50→40°C supply over T_ext -7→15°C. Off above 15°C, restarts at 13°C. ΔT = 5K constant.
# Cold: fixed 8°C supply, 13°C return. Off below 26°C, restarts at 27°C.
# Hysteresis prevents short-cycling.
# endregion

def T_hot_water_supply(T_ext: float, heating_active: bool) -> float | None:
    """Hot water supply temp [°C], or None if circuit off."""
    if heating_active:
        if T_ext >= T_HW_EXT_HIGH_C:
            return None
    else:
        if T_ext > T_HW_EXT_HYST_C:
            return None
    return float(np.interp(T_ext,
        [T_HW_EXT_LOW_C, T_HW_EXT_HIGH_C],
        [T_HW_SUPPLY_MAX, T_HW_SUPPLY_MIN]))


def T_cold_water_supply(T_ext: float, cooling_active: bool) -> float | None:
    """Cold water supply temp [°C], or None if circuit off. Fixed 8°C when on."""
    if cooling_active:
        if T_ext <= T_CW_EXT_LOW_C:
            return None
    else:
        if T_ext < T_CW_EXT_HYST_C:
            return None
    return T_CW_SUPPLY_C


# -----------------------------------------------------------------------------
# 5. POST-HOC HVAC DECOMPOSITION  (plotting only — never call from dT_dt)
# -----------------------------------------------------------------------------

def build_Q_hvac_array(
    T_in_array: np.ndarray,
    T_ext_array: np.ndarray,
    n_people_array: np.ndarray,
    Q_int_array: np.ndarray,
    dates: pd.DatetimeIndex,
) -> tuple:
    """
    Decompose HVAC into Q_heat, Q_cool, Q_vent + water flows.
    Returns: (Q_total, Q_heat, Q_cool, Q_vent,
              Q_water_heat, Q_water_cool, T_hw_arr, T_cw_arr,
              Q_air_m3s_arr, curtain_on)
    """
    n = len(T_in_array)
    Q_total      = np.zeros(n)
    Q_heat       = np.zeros(n)
    Q_cool       = np.zeros(n)
    Q_vent       = np.zeros(n)
    Q_water_heat = np.zeros(n)
    Q_water_cool = np.zeros(n)
    T_hw_arr     = np.full(n, np.nan)
    T_cw_arr     = np.full(n, np.nan)
    Q_air_m3s_arr = np.zeros(n)
    curtain_on   = np.zeros(n, dtype=bool)

    heating_active = False
    cooling_active = False

    for i, (T_in, T_ext, n_ppl, ts) in enumerate(
            zip(T_in_array, T_ext_array, n_people_array, dates)):

        hour      = ts.hour
        T_set     = T_setpoint(T_ext, hour)
        Q_air_m3s = airflow_total(n_ppl) / 3600.0
        T_mix     = FRAC_RETURN_AIR * T_in + (1.0 - FRAC_RETURN_AIR) * T_ext

        T_hw = T_hot_water_supply(T_ext, heating_active)
        T_cw = T_cold_water_supply(T_ext, cooling_active)
        heating_active = T_hw is not None
        cooling_active = T_cw is not None

        if T_hw is not None: T_hw_arr[i] = T_hw
        if T_cw is not None: T_cw_arr[i] = T_cw

        # Curtain state
        is_curtain = curtain_active(T_ext, hour)
        curtain_on[i] = is_curtain

        if np.isnan(T_set):
            # Dead band — pure ventilation with outdoor air
            q = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_ext) * (1 - FRAC_RETURN_AIR)
            Q_vent[i] = Q_total[i] = q

        elif T_set > T_in:
            # Heating — airflow boost if needed
            if T_in < T_set and heating_active:
                dtype  = _day_type(ts)
                is_night = 1 <= hour < 5
                offset = T_TUN_OFFSET_NIGHT_C if is_night else T_TUN_OFFSET_DAY_C
                T_tun  = T_ext + offset
                V_inf  = v_inf_m3s(hour, dtype)
                UA_eff = UA_FACADE_W_K + UA_TUN_WALL_W_K + RHO_CP_AIR_J_M3_K * V_inf
                Q_stair_flow = q_stair_m3s(T_ext, hour)

                Q_curtain_zone = 0.0
                if is_curtain:
                    Q_curtain_zone = ETA_IN_CURTAIN * RHO_CP_AIR_J_M3_K * Q_CURTAIN_AIR_M3S * DT_JET_K

                Q_gains = (UA_eff * (T_tun - T_in)
                           + UA_SOIL_W_K * (T_SOIL_C - T_in)
                           + Q_int_array[i]
                           + RHO_CP_AIR_J_M3_K * Q_stair_flow * (T_ext - T_in)
                           + Q_curtain_zone)

                heating_capacity = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_BLOW_HEAT_C - T_in)
                if Q_gains + heating_capacity < 0:
                    Q_air_required = -Q_gains / (RHO_CP_AIR_J_M3_K * (T_BLOW_HEAT_C - T_in))
                    Q_air_m3s = min(Q_air_required, AIRFLOW_MAX_M3H / 3600.0)

                q = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_BLOW_HEAT_C) # return air already accounted for in T_blow_heat_C
                Q_heat[i] = q
            else:
                q = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_ext) * (1 - FRAC_RETURN_AIR)
                Q_vent[i] = q
            Q_total[i] = q
            if heating_active and q != 0.0:
                Q_water_heat[i] = (Q_air_m3s * RHO_CP_AIR_J_M3_K * (T_BLOW_HEAT_C - T_mix)
                                   / (RHO_GLYCOL_KG_M3 * CP_GLYCOL_J_KG_K * DT_WATER_HEAT_K))

        else:
            # Cooling — airflow boost if hygiene flow can't meet load
            if cooling_active and T_in > T_BLOW_COOL_C:
                dtype  = _day_type(ts)
                is_night = 1 <= hour < 5
                offset = T_TUN_OFFSET_NIGHT_C if is_night else T_TUN_OFFSET_DAY_C
                T_tun  = T_ext + offset
                V_inf  = v_inf_m3s(hour, dtype)
                UA_eff = UA_FACADE_W_K + UA_TUN_WALL_W_K + RHO_CP_AIR_J_M3_K * V_inf
                Q_stair_flow = q_stair_m3s(T_ext, hour)

                Q_curtain_zone = 0.0
                if is_curtain:
                    Q_curtain_zone = ETA_IN_CURTAIN * RHO_CP_AIR_J_M3_K * Q_CURTAIN_AIR_M3S * DT_JET_K

                Q_load = (UA_eff * (T_tun - T_in)
                          + UA_SOIL_W_K * (T_SOIL_C - T_in)
                          + Q_int_array[i]
                          + RHO_CP_AIR_J_M3_K * Q_stair_flow * (T_ext - T_in)
                          + Q_curtain_zone)

                Q_cool_hygiene = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_BLOW_COOL_C)
                if Q_load > Q_cool_hygiene:
                    Q_air_required = Q_load / (RHO_CP_AIR_J_M3_K * (T_in - T_BLOW_COOL_C))
                    Q_air_m3s = min(Q_air_required, AIRFLOW_MAX_M3H / 3600.0)

                q = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_BLOW_COOL_C) # return air already accounted for in T_blow_cool_C
                Q_cool[i] = q
            else:
                q = RHO_CP_AIR_J_M3_K * Q_air_m3s * (T_in - T_ext) * (1 - FRAC_RETURN_AIR)
                Q_vent[i] = q
            Q_total[i] = q
            if cooling_active and q != 0.0:
                Q_water_cool[i] = (Q_air_m3s * RHO_CP_AIR_J_M3_K * (T_mix - T_BLOW_COOL_C)
                                   / (RHO_GLYCOL_KG_M3 * CP_GLYCOL_J_KG_K * DT_WATER_COOL_K))

        Q_air_m3s_arr[i] = Q_air_m3s

    return (Q_total, Q_heat, Q_cool, Q_vent,
            Q_water_heat, Q_water_cool, T_hw_arr, T_cw_arr,
            Q_air_m3s_arr, curtain_on)