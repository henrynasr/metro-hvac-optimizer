# =============================================================================
# humidity.py — Psychrometric utilities and indoor humidity computation
# Post-hoc only — never called from dT_dt.
# Functions: P_sat, W_from_T_RH, W_sat, RH_from_T_W, T_dew,
#            compute_humidity
# =============================================================================

import numpy as np
import pandas as pd

from constants import (
    RHO_AIR_KG_M3, FRAC_RETURN_AIR,
    T_CW_SUPPLY_C, T_CW_RETURN_C,
    WATTS_LATENT_PP,
    RH_TARGET_LOW, RH_TARGET_HIGH,
)
from occupancy import v_inf_m3s, _day_type
from regulation import q_stair_m3s, T_setpoint, T_cold_water_supply

# Standard atmospheric pressure [Pa]
P_ATM = 101_325.0

# Latent heat of vaporization [J/g]
H_FG_J_G = 2_450.0

# Coil surface temperature [°C] — average of cold water supply and return.
# Condensation happens on the fins at this temperature. [ASSUMPTION]
T_COIL_C = (T_CW_SUPPLY_C + T_CW_RETURN_C) / 2.0   # 10.5°C

# Latent emission per person [g/s]
M_DOT_PP_GS = WATTS_LATENT_PP / H_FG_J_G   # 60W / 2450 J/g ≈ 0.0245 g/s


# -----------------------------------------------------------------------------
# 1. PSYCHROMETRIC PRIMITIVES
# -----------------------------------------------------------------------------

def P_sat(T: float) -> float:
    """Saturation vapour pressure [Pa] from temperature [°C]. Magnus formula."""
    return 610.94 * np.exp(17.625 * T / (T + 243.04))


def W_from_T_RH(T: float, RH: float) -> float:
    """Humidity ratio [g/kg dry air] from temperature [°C] and RH [0–1]."""
    Ps = P_sat(T)
    Pw = RH * Ps
    return 622.0 * Pw / (P_ATM - Pw)


def W_sat(T: float) -> float:
    """Saturation humidity ratio [g/kg] at temperature T [°C]."""
    return W_from_T_RH(T, 1.0)


def RH_from_T_W(T: float, W: float) -> float:
    """Relative humidity [0–1] from temperature [°C] and humidity ratio [g/kg]."""
    Ps = P_sat(T)
    Pw = W * P_ATM / (622.0 + W)
    return np.clip(Pw / Ps, 0.0, 1.0)


def T_dew(W: float) -> float:
    """Dew point temperature [°C] from humidity ratio [g/kg]."""
    Pw = W * P_ATM / (622.0 + W)
    Pw = max(Pw, 1.0)   # guard against log(0)
    return 243.04 * np.log(Pw / 610.94) / (17.625 - np.log(Pw / 610.94))


# -----------------------------------------------------------------------------
# 2. INDOOR HUMIDITY COMPUTATION — post-hoc, one timestep
# -----------------------------------------------------------------------------

def _compute_W_in_step(
    T_ext: float,
    RH_ext: float,
    n_people: float,
    Q_air_m3s: float,
    Q_stair: float,
    V_inf: float,
    cooling_active: bool,
) -> float:
    """
    Solve moisture balance for W_in [g/kg] at one timestep.

    Heating / dead band / cooling without condensation:
        W_in = W_ext + M_people / (ρ × (Q_stair + V_inf + 0.30 × Q_air))

    Cooling with condensation (T_dew of mixed air > T_coil):
        W_supply = W_sat(T_coil), full balance solved for W_in.
    """
    W_ext = W_from_T_RH(T_ext, RH_ext)

    # Moisture from people [g/s]
    M_people = n_people * M_DOT_PP_GS

    # Fresh air flow — only the outdoor fraction of AHU dilutes moisture
    Q_fresh_ahu = (1.0 - FRAC_RETURN_AIR) * Q_air_m3s   # 0.30 × Q_air
    Q_fresh_total = Q_stair + V_inf + Q_fresh_ahu         # total outdoor air [m³/s]

    if Q_fresh_total < 1e-6:
        # No ventilation (shouldn't happen in practice) — can't solve, return outdoor
        return W_ext

    # Default: heating, dead band, or cooling without condensation
    W_in = W_ext + M_people / (RHO_AIR_KG_M3 * Q_fresh_total)

    # Check if cooling with condensation applies
    if cooling_active:
        # Mixed air humidity before coil
        W_mix = FRAC_RETURN_AIR * W_in + (1.0 - FRAC_RETURN_AIR) * W_ext
        T_dew_mix = T_dew(W_mix)

        if T_dew_mix > T_COIL_C:
            # Condensation on coil — W_supply is saturated at coil temp
            W_supply = W_sat(T_COIL_C)

            # Full balance: W_in × ρ × (Q_stair + V_inf + Q_air) =
            #   M_people + ρ × W_ext × (Q_stair + V_inf) + ρ × Q_air × W_supply
            Q_total = Q_stair + V_inf + Q_air_m3s
            numerator = (M_people
                         + RHO_AIR_KG_M3 * (Q_stair + V_inf) * W_ext
                         + RHO_AIR_KG_M3 * Q_air_m3s * W_supply)
            W_in = numerator / (RHO_AIR_KG_M3 * Q_total)

    return max(W_in, 0.0)


# -----------------------------------------------------------------------------
# 3. FULL HUMIDITY ARRAYS — loop over simulation results
# -----------------------------------------------------------------------------

def compute_humidity(r: dict) -> dict:
    """
    Compute indoor humidity from simulation results dict.

    Expects r to contain: T_in, T_ext, n_people, Q_air_m3s_arr, dates
    and requires RH_ext in the weather data.

    Returns dict with:
        W_ext, W_in, RH_in, W_supply arrays,
        comfort stats (% hours in/out of RH band).
    """
    dates     = r["dates"]
    T_in      = r["T_in"]
    T_ext     = r["T_ext"]
    RH_ext    = r["RH_ext"]
    n_people  = r["n_people"]
    Q_air_m3s = r["Q_air_m3s_arr"]

    n = len(dates)
    W_ext_arr    = np.empty(n)
    W_in_arr     = np.empty(n)
    RH_in_arr    = np.empty(n)
    W_supply_arr = np.empty(n)
    condensation = np.zeros(n, dtype=bool)
    latent_cool_W = np.zeros(n)

    cooling_active = False

    for i in range(n):
        ts    = dates[i]
        hour  = ts.hour
        dtype = _day_type(ts)

        # Flows at this timestep
        Q_stair = q_stair_m3s(T_ext[i], hour)
        V_inf   = v_inf_m3s(hour, dtype)

        # Water circuit state for cooling
        T_cw = T_cold_water_supply(T_ext[i], cooling_active)
        cooling_active = T_cw is not None

        # Outdoor humidity
        W_ext_arr[i] = W_from_T_RH(T_ext[i], RH_ext[i])

        # Solve for W_in
        W_in_arr[i] = _compute_W_in_step(
            T_ext[i], RH_ext[i], n_people[i],
            Q_air_m3s[i], Q_stair, V_inf, cooling_active,
        )

        # RH_in
        RH_in_arr[i] = RH_from_T_W(T_in[i], W_in_arr[i])

        # Condensation and latent load
        W_mix = (FRAC_RETURN_AIR * W_in_arr[i]
                 + (1.0 - FRAC_RETURN_AIR) * W_ext_arr[i])
        T_dew_mix = T_dew(W_mix)

        if cooling_active and T_dew_mix > T_COIL_C:
            condensation[i] = True
            W_supply_arr[i] = W_sat(T_COIL_C)
            # Latent load = mass of water removed × H_fg
            dW = W_mix - W_supply_arr[i]   # g/kg removed
            m_dot_condensate = RHO_AIR_KG_M3 * Q_air_m3s[i] * dW / 1000.0  # kg/s
            latent_cool_W[i] = m_dot_condensate * H_FG_J_G * 1000.0        # W
        else:
            W_supply_arr[i] = W_mix  # no condensation — air passes through unchanged

    # Comfort stats
    is_service = np.array([not (1 <= ts.hour < 5) for ts in dates])
    rh_ok      = (RH_in_arr >= RH_TARGET_LOW) & (RH_in_arr <= RH_TARGET_HIGH)
    rh_low     = RH_in_arr < RH_TARGET_LOW
    rh_high    = RH_in_arr > RH_TARGET_HIGH

    service_hours = is_service.sum()

    return {
        "W_ext":        W_ext_arr,
        "W_in":         W_in_arr,
        "RH_in":        RH_in_arr,
        "W_supply":     W_supply_arr,
        "condensation": condensation,
        "latent_cool_W": latent_cool_W,

        # Comfort summary
        "pct_RH_ok_service":   (rh_ok & is_service).sum() / service_hours * 100,
        "pct_RH_low_service":  (rh_low & is_service).sum() / service_hours * 100,
        "pct_RH_high_service": (rh_high & is_service).sum() / service_hours * 100,
        "hours_condensation":  condensation.sum(),
        "latent_cool_total_kWh": latent_cool_W.sum() / 1000.0,   # W × 1h = Wh, /1000 = kWh
    }

# -----------------------------------------------------------------------------
# 4. SMOKE TEST
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    # Quick sanity checks on psychrometric functions
    print("=== Psychrometric checks ===")
    print(f"P_sat(20°C)  = {P_sat(20):.0f} Pa   (expect ~2338)")
    print(f"W_sat(20°C)  = {W_sat(20):.2f} g/kg (expect ~14.7)")
    print(f"W(20°C, 50%) = {W_from_T_RH(20, 0.50):.2f} g/kg (expect ~7.3)")
    print(f"RH(20, 7.3)  = {RH_from_T_W(20, 7.3):.3f}       (expect ~0.50)")
    print(f"T_dew(7.3)   = {T_dew(7.3):.1f} °C    (expect ~9.2)")
    print(f"W_sat(T_coil={T_COIL_C}°C) = {W_sat(T_COIL_C):.2f} g/kg")
    print(f"M_dot/person = {M_DOT_PP_GS*3600:.1f} g/h  (expect ~88)")