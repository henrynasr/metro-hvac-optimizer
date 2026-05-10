# =============================================================================
# emissions.py — Electricity consumption and CO₂ emissions from HVAC
# Inputs: Q arrays from build_Q_hvac_array + RTE éCO2mix carbon intensity
# Outputs: P_elec_* [W], E_elec_* [kWh], CO2_* [kgCO₂], totals
# Post-hoc only — never called from dT_dt
# =============================================================================

import numpy as np
import pandas as pd

from constants import (
    COP_COOL, COP_HEAT,
    P_FAN_RATED_W, AIRFLOW_MAX_M3H,
    ELEC_PRICE_EUR_KWH,
)


# -----------------------------------------------------------------------------
# 1. RTE CO₂ INTENSITY  — load and align to simulation DatetimeIndex
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
# 2. ELECTRICAL POWER  — thermal load → electrical demand
# -----------------------------------------------------------------------------

def compute_elec_power(
    Q_heat: np.ndarray,
    Q_cool: np.ndarray,
    Q_vent: np.ndarray,
    Q_air_m3s_arr: np.ndarray,
) -> tuple:
    """
    Convert thermal powers [W] to electrical powers [W].

    Fan power uses the affinity cube law:
        P_fan = P_rated × (Q_actual / Q_rated)³
    This is the dominant savings lever — reducing flow by 20% cuts fan power ~49%.

    Returns: P_heat, P_cool, P_fan  [W electrical], each array len n
    """
    P_heat = np.abs(Q_heat) / COP_HEAT
    P_cool = np.abs(Q_cool) / COP_COOL

    # Fan runs whenever any HVAC mode is active (including dead band ventilation)
    fan_on = (np.abs(Q_heat) + np.abs(Q_cool) + np.abs(Q_vent)) > 0

    # Cube law: P_fan ∝ (Q/Q_max)³
    Q_rated_m3s = AIRFLOW_MAX_M3H / 3600.0
    flow_ratio  = Q_air_m3s_arr / Q_rated_m3s
    P_fan       = fan_on.astype(float) * P_FAN_RATED_W * flow_ratio**3

    return P_heat, P_cool, P_fan


# -----------------------------------------------------------------------------
# 3. ENERGY + CO₂  — integrate over time, multiply by carbon intensity
# -----------------------------------------------------------------------------

def compute_emissions(
    Q_heat: np.ndarray,
    Q_cool: np.ndarray,
    Q_vent: np.ndarray,
    Q_air_m3s_arr: np.ndarray,
    dates: pd.DatetimeIndex,
    co2_intensity: np.ndarray,
    dt_s: float = 3600.0,
) -> dict:
    """
    Full emissions pipeline. Returns dict with energy [kWh] and CO₂ [kgCO₂]
    arrays (per timestep) and annual totals.
    """
    P_heat, P_cool, P_fan = compute_elec_power(Q_heat, Q_cool, Q_vent, Q_air_m3s_arr)
    P_total = P_heat + P_cool + P_fan

    dt_h = dt_s / 3600.0

    E_heat  = P_heat  * dt_h / 1000.0
    E_cool  = P_cool  * dt_h / 1000.0
    E_fan   = P_fan   * dt_h / 1000.0
    E_total = P_total * dt_h / 1000.0

    CO2_heat  = E_heat  * co2_intensity / 1000.0
    CO2_cool  = E_cool  * co2_intensity / 1000.0
    CO2_fan   = E_fan   * co2_intensity / 1000.0
    CO2_total = E_total * co2_intensity / 1000.0

    return {
        "P_heat_W":    P_heat,
        "P_cool_W":    P_cool,
        "P_fan_W":     P_fan,
        "E_heat_kWh":  E_heat,
        "E_cool_kWh":  E_cool,
        "E_fan_kWh":   E_fan,
        "E_total_kWh": E_total,
        "CO2_heat_kg":  CO2_heat,
        "CO2_cool_kg":  CO2_cool,
        "CO2_fan_kg":   CO2_fan,
        "CO2_total_kg": CO2_total,
        "E_heat_total_kWh":   E_heat.sum(),
        "E_cool_total_kWh":   E_cool.sum(),
        "E_fan_total_kWh":    E_fan.sum(),
        "E_annual_kWh":       E_total.sum(),
        "CO2_annual_kgCO2":   CO2_total.sum(),
        "cost_annual_eur": E_total.sum() * ELEC_PRICE_EUR_KWH,
        "cost_heat_eur":   E_heat.sum()  * ELEC_PRICE_EUR_KWH,
        "cost_cool_eur":   E_cool.sum()  * ELEC_PRICE_EUR_KWH,
        "cost_fan_eur":    E_fan.sum()   * ELEC_PRICE_EUR_KWH,
    }