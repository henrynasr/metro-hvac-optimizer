# =============================================================================
# emissions.py — Electricity consumption and CO₂ emissions from HVAC
# Inputs: Q arrays from build_Q_hvac_array + RTE éCO2mix carbon intensity
# Outputs: P_elec_* [W], E_elec_* [kWh], CO2_* [kgCO₂], totals
# Post-hoc only — never called from dT_dt
# =============================================================================

import numpy as np
import pandas as pd

from constants import (
    COP_COOL, COP_HEAT, ETA_CARNOT_HEAT,
    P_FAN_RATED_W, AIRFLOW_MAX_M3H,
    ELEC_PRICE_EUR_KWH,
    RHO_CP_AIR_J_M3_K, Q_CURTAIN_AIR_M3S, DT_JET_K, P_FAN_CURTAIN_W,
)


# -----------------------------------------------------------------------------
# 1. HOROSAISONNIER TARIFF  — Tarif Jaune BT > 36 kVA, Courte Utilisation 2025
# -----------------------------------------------------------------------------
# Winter: Nov 1 – Mar 31.  Summer: Apr 1 – Oct 31.
# HC: 22h–06h Mon–Sat, all day Sunday.  HP: everything else.

TARIF_HPH = 0.215   # €/kWh HT — Heures Pleines Hiver
TARIF_HCH = 0.155   # €/kWh HT — Heures Creuses Hiver
TARIF_HPE = 0.118   # €/kWh HT — Heures Pleines Été
TARIF_HCE = 0.097   # €/kWh HT — Heures Creuses Été


def get_hourly_price(dates: pd.DatetimeIndex, flat: bool = False) -> np.ndarray:
    """
    Return €/kWh for each timestep.
    flat=True  → constant ELEC_PRICE_EUR_KWH (old behavior).
    flat=False → horosaisonnier 4-period tariff.
    """
    if flat:
        return np.full(len(dates), ELEC_PRICE_EUR_KWH)

    month = dates.month
    hour  = dates.hour
    dow   = dates.dayofweek          # 0=Mon, 6=Sun

    is_winter = (month >= 11) | (month <= 3)
    is_hc     = (dow == 6) | (hour >= 22) | (hour < 6)

    return np.where(
        is_winter,
        np.where(is_hc, TARIF_HCH, TARIF_HPH),
        np.where(is_hc, TARIF_HCE, TARIF_HPE),
    )


# -----------------------------------------------------------------------------
# 2. RTE CO₂ INTENSITY  — load and align to simulation DatetimeIndex
# -----------------------------------------------------------------------------

def load_co2_intensity(rte_path: str, dates: pd.DatetimeIndex) -> np.ndarray:
    """
    Parse RTE éCO2mix file, return carbon intensity [gCO₂/kWh] aligned to dates.
    RTE file is 30-min resolution with alternating empty rows — forward-filled to hourly.
    """
    df = pd.read_csv(rte_path, sep="\t", encoding="ISO-8859-1",
                     usecols=["Date", "Heures", "Taux de Co2"])

    df["datetime"] = pd.to_datetime(
        df["Date"] + " " + df["Heures"], dayfirst=False
    )

    df = df.dropna(subset=["Taux de Co2"]).set_index("datetime")
    df = df["Taux de Co2"].sort_index()

    df_aligned = df.reindex(dates, method="nearest", tolerance="1h")
    df_aligned = df_aligned.ffill().bfill()

    return df_aligned.values.astype(float)


# -----------------------------------------------------------------------------
# 3. ELECTRICAL POWER  — thermal load → electrical demand
# -----------------------------------------------------------------------------

def compute_elec_power(
    Q_heat: np.ndarray,
    Q_cool: np.ndarray,
    Q_vent: np.ndarray,
    Q_air_m3s_arr: np.ndarray,
    curtain_on: np.ndarray,
    T_ext_array: np.ndarray,
    T_hw_arr: np.ndarray,
) -> tuple:
    """
    Convert thermal powers [W] to electrical powers [W].

    Heating COP is variable — Carnot-based, air-source heat pump:
        COP(T_hw, T_ext) = η_carnot × T_hw_K / (T_hw_K − T_ext_K)
    Falls back to COP_HEAT constant when hot water circuit is off (T_hw = NaN).
    Same variable COP applies to air curtain (parallel branch on same HP).

    Fan power uses the affinity cube law:
        P_fan = P_rated × (Q_actual / Q_rated)³

    Returns: P_heat, P_cool, P_fan, P_curtain  [W electrical], each array len n
    """
    # Variable COP for heating — air-source HP
    T_hw_K  = T_hw_arr + 273.15
    T_ext_K = T_ext_array + 273.15
    dT_K    = T_hw_K - T_ext_K

    # Guard against division by zero (T_hw ≈ T_ext) and NaN (circuit off)
    safe_dT = np.where((dT_K > 1.0) & np.isfinite(T_hw_K), dT_K, np.nan)
    COP_heat_var = ETA_CARNOT_HEAT * T_hw_K / safe_dT

    # Clamp to physical range [2.0, 8.0] and fallback to constant when NaN
    COP_heat_var = np.clip(COP_heat_var, 2.0, 8.0)
    COP_heat_var = np.where(np.isfinite(COP_heat_var), COP_heat_var, COP_HEAT)

    P_heat = np.abs(Q_heat) / COP_heat_var
    P_cool = np.abs(Q_cool) / COP_COOL

    # Fan runs whenever any HVAC mode is active (including dead band ventilation)
    fan_on = (np.abs(Q_heat) + np.abs(Q_cool) + np.abs(Q_vent)) > 0

    # Cube law: P_fan ∝ (Q/Q_max)³
    Q_rated_m3s = AIRFLOW_MAX_M3H / 3600.0
    flow_ratio  = Q_air_m3s_arr / Q_rated_m3s
    P_fan       = fan_on.astype(float) * P_FAN_RATED_W * flow_ratio**3

    # Air curtain — same HP, same variable COP
    P_heat_curtain = RHO_CP_AIR_J_M3_K * Q_CURTAIN_AIR_M3S * DT_JET_K  # 8442 W thermal
    P_curtain = curtain_on.astype(float) * (P_heat_curtain / COP_heat_var + P_FAN_CURTAIN_W)

    return P_heat, P_cool, P_fan, P_curtain


# -----------------------------------------------------------------------------
# 4. ENERGY + CO₂ + COST  — integrate over time, apply tariff
# -----------------------------------------------------------------------------

def compute_emissions(
    Q_heat: np.ndarray,
    Q_cool: np.ndarray,
    Q_vent: np.ndarray,
    Q_air_m3s_arr: np.ndarray,
    curtain_on: np.ndarray,
    co2_intensity: np.ndarray,
    T_ext_array: np.ndarray,
    T_hw_arr: np.ndarray,
    dates: pd.DatetimeIndex = None,
    flat_tariff: bool = False,
    dt_s: float = 3600.0,
) -> dict:
    """
    Full emissions pipeline. Returns dict with energy [kWh], CO₂ [kgCO₂],
    and cost [€] arrays (per timestep) and annual totals.

    Tariff: horosaisonnier by default (4-period: HPH/HCH/HPE/HCE).
    Pass flat_tariff=True for constant ELEC_PRICE_EUR_KWH (old behavior).
    """
    P_heat, P_cool, P_fan, P_curtain = compute_elec_power(
        Q_heat, Q_cool, Q_vent, Q_air_m3s_arr, curtain_on,
        T_ext_array, T_hw_arr)
    P_total = P_heat + P_cool + P_fan + P_curtain

    dt_h = dt_s / 3600.0

    E_heat    = P_heat    * dt_h / 1000.0
    E_cool    = P_cool    * dt_h / 1000.0
    E_fan     = P_fan     * dt_h / 1000.0
    E_curtain = P_curtain * dt_h / 1000.0
    E_total   = P_total   * dt_h / 1000.0

    CO2_heat    = E_heat    * co2_intensity / 1000.0
    CO2_cool    = E_cool    * co2_intensity / 1000.0
    CO2_fan     = E_fan     * co2_intensity / 1000.0
    CO2_curtain = E_curtain * co2_intensity / 1000.0
    CO2_total   = E_total   * co2_intensity / 1000.0

    # --- Hourly tariff ---
    if dates is not None:
        price = get_hourly_price(dates, flat=flat_tariff)
    else:
        price = np.full(len(E_total), ELEC_PRICE_EUR_KWH)

    cost_heat    = E_heat    * price
    cost_cool    = E_cool    * price
    cost_fan     = E_fan     * price
    cost_curtain = E_curtain * price
    cost_total   = E_total   * price

    return {
        "P_heat_W":    P_heat,
        "P_cool_W":    P_cool,
        "P_fan_W":     P_fan,
        "P_curtain_W": P_curtain,
        "E_heat_kWh":    E_heat,
        "E_cool_kWh":    E_cool,
        "E_fan_kWh":     E_fan,
        "E_curtain_kWh": E_curtain,
        "E_total_kWh":   E_total,
        "CO2_heat_kg":    CO2_heat,
        "CO2_cool_kg":    CO2_cool,
        "CO2_fan_kg":     CO2_fan,
        "CO2_curtain_kg": CO2_curtain,
        "CO2_total_kg":   CO2_total,
        "price_eur_kwh":   price,
        "cost_heat_eur_ts":    cost_heat,
        "cost_cool_eur_ts":    cost_cool,
        "cost_fan_eur_ts":     cost_fan,
        "cost_curtain_eur_ts": cost_curtain,
        "cost_total_eur_ts":   cost_total,
        "E_heat_total_kWh":    E_heat.sum(),
        "E_cool_total_kWh":    E_cool.sum(),
        "E_fan_total_kWh":     E_fan.sum(),
        "E_curtain_total_kWh": E_curtain.sum(),
        "E_annual_kWh":        E_total.sum(),
        "CO2_annual_kgCO2":    CO2_total.sum(),
        "cost_annual_eur":   cost_total.sum(),
        "cost_heat_eur":     cost_heat.sum(),
        "cost_cool_eur":     cost_cool.sum(),
        "cost_fan_eur":      cost_fan.sum(),
        "cost_curtain_eur":  cost_curtain.sum(),
    }