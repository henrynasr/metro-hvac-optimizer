# =============================================================================
# thermal_model.py — Single-run thermal simulation, one platform zone
# =============================================================================

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from constants import (
    C_TOTAL_J_K, T_SOIL_C,
    T_TUN_OFFSET_DAY_C, T_TUN_OFFSET_NIGHT_C,
    UA_FACADE_W_K, UA_TUN_WALL_W_K, UA_SOIL_W_K, RHO_CP_AIR_J_M3_K,
    Q_STAIR_M3S, AIRFLOW_MAX_M3H
)
from occupancy import load_profiles, build_Q_array
from regulation import T_setpoint, airflow_total, dT_dt, build_Q_hvac_array
from utils import load_data, style_axes
from emissions import load_co2_intensity, compute_emissions

# -----------------------------------------------------------------------------
# 1. INPUTS
# -----------------------------------------------------------------------------

SIM_START   = "2024-01-01"
SIM_END     = "2024-12-31"
XLSX_PATH   = "data/raw/Defense_Occupation_Normalised.xlsx"
WEATHER_CSV = "data/raw/paris_weather.csv"
T0          = 25.0   # °C — initial T_in

df          = load_data(WEATHER_CSV)
df_sim      = df.loc[SIM_START:SIM_END].copy()
dates       = df_sim.index
T_ext_array = df_sim["temperature_2m"].values
t_array     = np.arange(len(dates)) * 3600.0

profiles             = load_profiles(XLSX_PATH)
Q_int_array, n_array = build_Q_array(dates, profiles)

# -----------------------------------------------------------------------------
# 2. ODE SOLVE
# -----------------------------------------------------------------------------

water_state = {"heating": False, "cooling": False}

sol = solve_ivp(
    fun=dT_dt,
    t_span=(t_array[0], t_array[-1]),
    y0=[T0],
    t_eval=t_array,
    args=(t_array, T_ext_array, Q_int_array, n_array, dates, water_state),
    method="RK45",
    max_step=3600.0,
)

T_in_array  = sol.y[0]
t_hours     = t_array / 3600.0

# T_tun per hour (day/night split for plot)
T_tun_array = np.array([
    T_ext_array[i] + (T_TUN_OFFSET_DAY_C if 5 <= dates[i].hour <= 23 else T_TUN_OFFSET_NIGHT_C)
    for i in range(len(dates))
])

Q_total, Q_heat, Q_cool, Q_vent, Q_water_heat, Q_water_cool, T_hw_arr, T_cw_arr = build_Q_hvac_array(
    T_in_array, T_ext_array, n_array, dates
)

RTE_PATH = "data/raw/eco2mix_2024.csv"

Q_air_m3s = AIRFLOW_MAX_M3H / 3600.0  # design flow for fan power calc

co2_intensity = load_co2_intensity(RTE_PATH, dates)
em = compute_emissions(Q_heat, Q_cool, Q_vent, Q_air_m3s, dates, co2_intensity)

# -----------------------------------------------------------------------------
# 3. PLOT (4×2 panel)
# -----------------------------------------------------------------------------

fig, axes = plt.subplots(4, 2, figsize=(16, 14), sharex=True)
fig.suptitle(
    f"Thermal model — Platform zone, full year 2024\n"
    f"facade→T_tun | concrete→T_soil | stair infiltration | AHU regulation",
    fontsize=13, y=0.98
)

# Left column
ax = axes[0, 0]
ax.plot(t_hours, T_ext_array, label="T_ext (outdoor)",  color="steelblue",  lw=1.2)
ax.plot(t_hours, T_tun_array, label="T_tun (tunnel)",   color="darkorange", lw=1.2, ls="--")
ax.axhline(T_SOIL_C, color="sienna", lw=1.0, ls=":", label=f"T_soil = {T_SOIL_C}°C")
ax.plot(t_hours, T_in_array,  label="T_in (platform)",  color="crimson",    lw=1.6)
style_axes(ax, title="Temperatures", ylabel="°C")
ax.legend(fontsize=9)

ax = axes[1, 0]
ax.plot(t_hours, Q_int_array / 1000, color="darkorchid", lw=1.4)
style_axes(ax, title="Internal heat gain (occupancy + equipment)", ylabel="kW")

ax = axes[2, 0]
ax.plot(t_hours, n_array, color="teal", lw=1.4)
style_axes(ax, title="Headcount", ylabel="persons")

ax = axes[3, 0]
ax.plot(t_hours, Q_total / 1000, color="black", lw=1.6, label="Q_hvac total")
ax.axhline(0, color="gray", lw=0.8, ls="--")
style_axes(ax, title="Total HVAC power (+ = cooling load on zone)", ylabel="kW", xlabel="Time [h]")
ax.legend(fontsize=9)

# Right column
ax = axes[0, 1]
ax.plot(t_hours, Q_heat / 1000, color="tomato", lw=1.4, label="Q_heat")
style_axes(ax, title="Heating power (negative = adds heat)", ylabel="kW")
ax.legend(fontsize=9)

ax = axes[1, 1]
ax.plot(t_hours, Q_cool / 1000, color="dodgerblue", lw=1.4, label="Q_cool")
style_axes(ax, title="Cooling power", ylabel="kW")
ax.legend(fontsize=9)

ax = axes[2, 1]
ax.plot(t_hours, Q_vent / 1000, color="seagreen", lw=1.4, label="Q_vent (dead band)")
ax.axhline(0, color="gray", lw=0.8, ls="--")
style_axes(ax, title="Ventilation heat exchange (dead band)", ylabel="kW")
ax.legend(fontsize=9)

# Steady-state annotation
UA_static   = UA_FACADE_W_K + UA_TUN_WALL_W_K + UA_SOIL_W_K
T_tun_mean  = T_ext_array.mean() + T_TUN_OFFSET_DAY_C
T_eq_approx = (
    (UA_FACADE_W_K + UA_TUN_WALL_W_K) * T_tun_mean
    + UA_SOIL_W_K * T_SOIL_C
    + Q_int_array.mean()
    + RHO_CP_AIR_J_M3_K * Q_STAIR_M3S * T_ext_array.mean()
) / (UA_static + RHO_CP_AIR_J_M3_K * Q_STAIR_M3S)

lines = [
    "Steady-state check (no HVAC, mean conditions):",
    f"  UA_facade+wall = {UA_FACADE_W_K + UA_TUN_WALL_W_K:.0f} W/K",
    f"  UA_soil        = {UA_SOIL_W_K:.0f} W/K",
    f"  T_tun_mean     = {T_tun_mean:.1f} °C",
    f"  T_soil         = {T_SOIL_C:.1f} °C",
    f"  Q_int_mean     = {Q_int_array.mean()/1000:.1f} kW",
    f"  → T_eq ≈ {T_eq_approx:.1f} °C (no HVAC)",
    "",
    f"  C_total = {C_TOTAL_J_K:.2e} J/K",
    f"  τ = C / UA = {C_TOTAL_J_K / UA_static / 3600:.0f} h",
]

axes[3, 1].set_visible(False)
ax_note = fig.add_axes([0.53, 0.02, 0.42, 0.22])
ax_note.axis("off")
ax_note.text(0.05, 0.95, "\n".join(lines), transform=ax_note.transAxes,
             va="top", ha="left", fontsize=17, fontfamily="monospace",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("images/thermal_model_2024.png", dpi=150)
plt.close()

print(f"Saved: images/thermal_model.png")
print(f"Simulation: {SIM_START} → {SIM_END}")
print(f"T_in:  {T_in_array.min():.1f}–{T_in_array.max():.1f} °C")
print(f"T_ext: {T_ext_array.min():.1f}–{T_ext_array.max():.1f} °C")
print(f"Q_int: {Q_int_array.min()/1e3:.1f}–{Q_int_array.max()/1e3:.1f} kW")
print(f"Q_hvac: {Q_total.min()/1e3:.1f}–{Q_total.max()/1e3:.1f} kW")
print(f"Steady-state T_eq (no HVAC): {T_eq_approx:.1f} °C")

print(f"\n=== Annual Energy & Emissions ===")
print(f"E_heating : {em['E_heat_total_kWh']:>8.0f} kWh")
print(f"E_cooling : {em['E_cool_total_kWh']:>8.0f} kWh")
print(f"E_fans    : {em['E_fan_total_kWh']:>8.0f} kWh")
print(f"E_TOTAL   : {em['E_annual_kWh']:>8.0f} kWh")
print(f"CO2_TOTAL : {em['CO2_annual_kgCO2']:>8.0f} kgCO₂")

# Monthly energy breakdown
monthly = pd.DataFrame({
    "E_heat_kWh":  em["E_heat_kWh"],
    "E_cool_kWh":  em["E_cool_kWh"],
    "E_fan_kWh":   em["E_fan_kWh"],
    "CO2_total_kg": em["CO2_total_kg"],
    "intensity":   co2_intensity,
}, index=dates)

monthly_sum = monthly.resample("ME").sum()
monthly_sum["intensity_mean"] = monthly.resample("ME")["intensity"].mean()

print("\n=== Monthly breakdown ===")
print(monthly_sum[["E_heat_kWh","E_cool_kWh","E_fan_kWh","CO2_total_kg", "intensity_mean"]].round(0).to_string())

print(f"\n=== Annual Cost ===")
print(f"Heating : {em['cost_heat_eur']:>8.0f} €")
print(f"Cooling : {em['cost_cool_eur']:>8.0f} €")
print(f"Fans    : {em['cost_fan_eur']:>8.0f} €")
print(f"TOTAL   : {em['cost_annual_eur']:>8.0f} €")