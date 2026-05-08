# =============================================================================
# thermal_model.py — Single-run thermal simulation for one platform zone
#
# Drives the lumped-capacitance ODE with:
#   - Real Paris weather (Open-Meteo ERA5)
#   - Real RATP occupancy profiles (Pôle La Défense)
#   - Full regulation layer (setpoint law, airflow modulation, HVAC power)
#   - Split envelope: facade (→ T_tun) + structural (→ T_soil)
#   - Dynamic infiltration: V̇_inf(hour, day_type) from occupancy.py
#
# Output: 4×2 panel plot saved to images/thermal_model.png
# =============================================================================

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from constants import (
    C_TOTAL_J_K,
    T_SOIL_C,
    T_TUN_OFFSET_C,
    UA_FACADE_W_K,
    UA_SOIL_W_K,
    RHO_CP_AIR_J_M3_K,
)
from occupancy import load_profiles, build_Q_array
from regulation import T_setpoint, airflow_total, dT_dt, build_Q_hvac_array
from utils import load_data, style_axes

# -----------------------------------------------------------------------------
# SECTION 1 — Inputs
# -----------------------------------------------------------------------------

SIM_START   = "2024-07-01"
SIM_END     = "2024-07-07"
XLSX_PATH   = "data/raw/Defense_Occupation_Normalised.xlsx"
WEATHER_CSV = "data/raw/paris_weather.csv"
T0          = 23.0   # °C — initial indoor temperature

df      = load_data(WEATHER_CSV)
df_sim  = df.loc[SIM_START:SIM_END].copy()
dates   = df_sim.index

T_ext_array = df_sim["temperature_2m"].values
t_array     = np.arange(len(dates)) * 3600.0   # seconds from t=0

profiles              = load_profiles(XLSX_PATH)
Q_int_array, n_array  = build_Q_array(dates, profiles)

# -----------------------------------------------------------------------------
# SECTION 2 — ODE solve
# -----------------------------------------------------------------------------

sol = solve_ivp(
    fun=dT_dt,
    t_span=(t_array[0], t_array[-1]),
    y0=[T0],
    t_eval=t_array,
    args=(t_array, T_ext_array, Q_int_array, n_array, dates),
    method="RK45",
    max_step=3600.0,
)

T_in_array = sol.y[0]

# Derived boundary arrays for plotting
T_tun_array  = T_ext_array + T_TUN_OFFSET_C

# Post-hoc HVAC decomposition
Q_total, Q_heat, Q_cool, Q_vent, Q_water_heat, Q_water_cool, T_hw_arr, T_cw_arr = build_Q_hvac_array(
    T_in_array, T_ext_array, n_array, dates
)

# Time axis in hours for plotting
t_hours = t_array / 3600.0

# -----------------------------------------------------------------------------
# SECTION 3 — Plot (4×2 panel)
# -----------------------------------------------------------------------------

fig, axes = plt.subplots(4, 2, figsize=(16, 14), sharex=True)
fig.suptitle(
    f"Thermal model — Platform zone, {SIM_START} to {SIM_END}\n"
    f"One buried zone | facade→T_tun | concrete→T_soil | infiltration V̇_inf(hour)",
    fontsize=13, y=0.98
)

# --- Left column: primary state variables ---

ax = axes[0, 0]
ax.plot(t_hours, T_ext_array,  label="T_ext (outdoor)",  color="steelblue",  lw=1.2)
ax.plot(t_hours, T_tun_array,  label="T_tun (tunnel)",   color="darkorange",  lw=1.2, ls="--")
ax.axhline(T_SOIL_C, color="sienna", lw=1.0, ls=":", label=f"T_soil = {T_SOIL_C}°C")
ax.plot(t_hours, T_in_array,   label="T_in (platform)",  color="crimson",     lw=1.6)
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
style_axes(ax, title="Total HVAC power (+ = cooling load on zone)",
           ylabel="kW", xlabel="Time [h]")
ax.axhline(0, color="gray", lw=0.8, ls="--")
ax.legend(fontsize=9)

# --- Right column: HVAC decomposition ---

ax = axes[0, 1]
ax.plot(t_hours, Q_heat / 1000, color="tomato",  lw=1.4, label="Q_heat (active heating)")
style_axes(ax, title="Heating power (negative = adds heat)", ylabel="kW")
ax.legend(fontsize=9)

ax = axes[1, 1]
ax.plot(t_hours, Q_cool / 1000, color="dodgerblue", lw=1.4, label="Q_cool (active cooling)")
style_axes(ax, title="Cooling power", ylabel="kW")
ax.legend(fontsize=9)

ax = axes[2, 1]
ax.plot(t_hours, Q_vent / 1000, color="seagreen", lw=1.4, label="Q_vent (dead-band ventilation)")
style_axes(ax, title="Ventilation heat exchange (dead band)", ylabel="kW")
ax.axhline(0, color="gray", lw=0.8, ls="--")
ax.legend(fontsize=9)

# Bottom-right: steady-state sanity check annotation
ax = axes[3, 1]
ax.axis("off")
UA_static = UA_FACADE_W_K + UA_SOIL_W_K
T_eq_approx = (
    UA_FACADE_W_K * (T_ext_array.mean() + T_TUN_OFFSET_C)
    + UA_SOIL_W_K * T_SOIL_C
    + Q_int_array.mean()
) / UA_static

lines = [
    "Steady-state check (no HVAC, mean conditions):",
    f"  UA_facade = {UA_FACADE_W_K:.0f} W/K",
    f"  UA_soil   = {UA_SOIL_W_K:.0f} W/K",
    f"  T_tun_mean = {T_ext_array.mean() + T_TUN_OFFSET_C:.1f} °C",
    f"  T_soil    = {T_SOIL_C:.1f} °C",
    f"  Q_int_mean = {Q_int_array.mean()/1000:.1f} kW",
    f"  → T_eq ≈ {T_eq_approx:.1f} °C (no cooling)",
    "",
    f"  C_total   = {C_TOTAL_J_K:.2e} J/K",
    f"  τ_static  = C / UA_static = {C_TOTAL_J_K / UA_static / 3600:.0f} h",
]
ax.text(0.05, 0.95, "\n".join(lines), transform=ax.transAxes,
        va="top", ha="left", fontsize=9,
        fontfamily="monospace",
        bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

plt.tight_layout(rect=[0, 0, 1, 0.96])
axes[3, 1].set_visible(False)
ax_note = fig.add_axes([0.53, 0.02, 0.42, 0.22])  # [left, bottom, width, height]
ax_note.axis("off")
ax_note.text(0.05, 0.95, "\n".join(lines), transform=ax_note.transAxes,
             va="top", ha="left", fontsize=17, fontfamily="monospace",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))
plt.savefig("images/thermal_model.png", dpi=150)
plt.close()
print("Saved: images/thermal_model.png")

# Quick console summary
print(f"\nSimulation: {SIM_START} → {SIM_END}")
print(f"T_in: {T_in_array.min():.1f}–{T_in_array.max():.1f} °C")
print(f"T_ext: {T_ext_array.min():.1f}–{T_ext_array.max():.1f} °C")
print(f"T_tun: {T_tun_array.min():.1f}–{T_tun_array.max():.1f} °C")
print(f"Q_int: {Q_int_array.min()/1e3:.1f}–{Q_int_array.max()/1e3:.1f} kW")
print(f"Q_hvac total: {Q_total.min()/1e3:.1f}–{Q_total.max()/1e3:.1f} kW")
print(f"Steady-state T_eq (no HVAC): {T_eq_approx:.1f} °C")