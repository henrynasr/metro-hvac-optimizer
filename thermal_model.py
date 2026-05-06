"""
thermal_model.py — single-run lumped-capacitance thermal model
for one metro station zone, driven by Paris weather and RATP
occupancy. Produces a three-panel diagnostic plot.
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

from utils import load_data, style_axes
from occupancy import load_profiles, build_Q_array
from regulation import dT_dt, build_Q_hvac_array


# ---------------------------------------------------------------------
# 1. Inputs — weather and occupancy over the simulation window
# ---------------------------------------------------------------------

df = load_data("data/raw/paris_weather.csv").sort_index()
df_sliced = df.loc["2024-08-05":"2024-08-12"]

t_array     = np.arange(len(df_sliced)) * 3600       # seconds
T_ext_array = df_sliced["temperature_2m"].values     # °C
dates       = df_sliced.index

profiles = load_profiles("data/raw/Defense_Occupation_Normalised.xlsx")
Q_array, N_people_array = build_Q_array(dates, profiles)

# ---------------------------------------------------------------------
# 2. Physics parameters, ODE solve and HVAC contribution
# ---------------------------------------------------------------------

UA = 5e3       # envelope conductance (W/K)
C  = 5e7       # lumped thermal capacitance (J/K)
T0 = [16.0]    # initial indoor temperature (°C)

sol = solve_ivp(
    fun=dT_dt,
    t_span=(t_array[0], t_array[-1]),
    y0=T0,
    t_eval=t_array,
    args=(t_array, T_ext_array, Q_array, N_people_array, UA, C),
)

Q_hvac_array, Q_heat_array, Q_cool_array, Q_vent_array = build_Q_hvac_array(sol.y[0], T_ext_array, N_people_array)


# ---------------------------------------------------------------------
# 3. Three-panel diagnostic plot
#    Top    : T_ext vs T_in
#    Middle : internal heat load Q
#    Bottom : occupant headcount n
# ---------------------------------------------------------------------

if __name__ == "__main__":

    hours = t_array / 3600
    fig, axes = plt.subplots(4, 2, figsize=(11, 10), sharex=True)
    fig.suptitle("Paris station thermal model — Second week of Aug 2024", fontsize = 13)
    # Top panel — outdoor vs indoor temperature
    axes[0,0].plot(hours, T_ext_array, label="T_ext (outdoor)")
    axes[0,0].plot(hours, sol.y[0],    label="T_in (indoor)")
    style_axes(axes[0,0], "", "", "Temperature (°C)")
    axes[0,0].legend(fontsize=12)

    # Second panel — internal heat load
    axes[1,0].plot(hours, Q_array / 1000, color="tab:orange")
    style_axes(axes[1,0], "", "", "Q internal (kW)")

    # Third panel — headcount
    axes[2,0].plot(hours, N_people_array, color="tab:green")
    style_axes(axes[2,0], "", "", "Headcount (persons)")

    # Bottom panel — HVAC contribution
    axes[3,0].plot(hours, Q_hvac_array / 1000, color="tab:blue")
    style_axes(axes[3,0], "", "Hours since start", "Q_HVAC (KW)")

    axes[0,1].plot(t_array / 3600, Q_heat_array / 1000, color="red")
    style_axes(axes[0,1], ylabel="Q_heat (kW)")

    axes[1,1].plot(t_array / 3600, Q_cool_array / 1000, color="blue")
    style_axes(axes[1,1], ylabel="Q_cool (kW)")

    axes[2,1].plot(t_array / 3600, Q_vent_array / 1000, color="grey")
    style_axes(axes[2,1], ylabel="Q_vent (kW)")

    plt.tight_layout()
    plt.savefig("images/thermal_model_aug2024.png", dpi=150)
    plt.close()