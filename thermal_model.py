# =============================================================================
# thermal_model.py — Single-run thermal simulation, one platform zone
# run_simulation() returns results dict — reusable by sweep/optimizer.
# plot_results() saves 4×2 panel plot.
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
    Q_STAIR_M3S,
)
from occupancy import load_profiles, build_Q_array
from regulation import dT_dt, build_Q_hvac_array
from utils import load_data, style_axes
from emissions import load_co2_intensity, compute_emissions


# -----------------------------------------------------------------------------
# 1. SIMULATION ENGINE
# -----------------------------------------------------------------------------

def run_simulation(
    sim_start: str = "2024-01-01",
    sim_end: str = "2024-12-31",
    T0: float = 25.0,
    weather_csv: str = "data/raw/paris_weather.csv",
    xlsx_path: str = "data/raw/Defense_Occupation_Normalised.xlsx",
    rte_path: str = "data/raw/eco2mix_2024.csv",
) -> dict:
    """
    Run full-year ODE + post-hoc HVAC decomposition + emissions.
    Returns dict with all arrays and summary scalars.
    """
    df       = load_data(weather_csv)
    df_sim   = df.loc[sim_start:sim_end].copy()
    dates    = df_sim.index
    T_ext    = df_sim["temperature_2m"].values
    t_array  = np.arange(len(dates)) * 3600.0

    profiles        = load_profiles(xlsx_path)
    Q_int, n_people = build_Q_array(dates, profiles)

    # ODE solve
    water_state = {"heating": False, "cooling": False}
    sol = solve_ivp(
        fun=dT_dt,
        t_span=(t_array[0], t_array[-1]),
        y0=[T0],
        t_eval=t_array,
        args=(t_array, T_ext, Q_int, n_people, dates, water_state),
        method="RK45",
        max_step=3600.0,
    )
    T_in = sol.y[0]

    # T_tun per hour — night = 01h–05h
    T_tun = np.array([
        T_ext[i] + (T_TUN_OFFSET_NIGHT_C if 1 <= dates[i].hour < 5 else T_TUN_OFFSET_DAY_C)
        for i in range(len(dates))
    ])

    # Post-hoc HVAC decomposition
    (Q_total, Q_heat, Q_cool, Q_vent,
     Q_water_heat, Q_water_cool, T_hw, T_cw,
     Q_air_m3s_arr, curtain_on) = build_Q_hvac_array(T_in, T_ext, n_people, Q_int, dates)

    # Emissions
    co2_intensity = load_co2_intensity(rte_path, dates)
    em = compute_emissions(Q_heat, Q_cool, Q_vent, Q_air_m3s_arr, curtain_on, dates, co2_intensity)

    return {
        "dates": dates, "t_array": t_array,
        "T_ext": T_ext, "T_in": T_in, "T_tun": T_tun,
        "Q_int": Q_int, "n_people": n_people,
        "Q_total": Q_total, "Q_heat": Q_heat, "Q_cool": Q_cool, "Q_vent": Q_vent,
        "Q_water_heat": Q_water_heat, "Q_water_cool": Q_water_cool,
        "T_hw": T_hw, "T_cw": T_cw,
        "Q_air_m3s_arr": Q_air_m3s_arr,
        "curtain_on": curtain_on,
        "co2_intensity": co2_intensity,
        "em": em,
    }


# -----------------------------------------------------------------------------
# 2. PLOT
# -----------------------------------------------------------------------------

def plot_results(r: dict, filename: str = "images/thermal_model_2024.png"):
    """4×2 panel plot from simulation results dict."""
    t_h = r["t_array"] / 3600.0

    fig, axes = plt.subplots(4, 2, figsize=(16, 14), sharex=True)
    fig.suptitle(
        "Thermal model — Platform zone, full year 2024\n"
        "facade→T_tun | concrete→T_soil | stair infiltration | AHU regulation",
        fontsize=13, y=0.98,
    )

    # Left column
    ax = axes[0, 0]
    ax.plot(t_h, r["T_ext"], label="T_ext",  color="steelblue",  lw=1.2)
    ax.plot(t_h, r["T_tun"], label="T_tun",  color="darkorange", lw=1.2, ls="--")
    ax.axhline(T_SOIL_C, color="sienna", lw=1.0, ls=":", label=f"T_soil = {T_SOIL_C}°C")
    ax.plot(t_h, r["T_in"],  label="T_in",   color="crimson",    lw=1.6)
    style_axes(ax, title="Temperatures", ylabel="°C")
    ax.legend(fontsize=9)

    ax = axes[1, 0]
    ax.plot(t_h, r["Q_int"] / 1000, color="darkorchid", lw=1.4)
    style_axes(ax, title="Internal heat gain", ylabel="kW")

    ax = axes[2, 0]
    ax.plot(t_h, r["n_people"], color="teal", lw=1.4)
    style_axes(ax, title="Headcount", ylabel="persons")

    ax = axes[3, 0]
    ax.plot(t_h, r["Q_total"] / 1000, color="black", lw=1.6, label="Q_hvac total")
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    style_axes(ax, title="Total HVAC power (+ = cooling)", ylabel="kW", xlabel="Time [h]")
    ax.legend(fontsize=9)

    # Right column
    ax = axes[0, 1]
    ax.plot(t_h, r["Q_heat"] / 1000, color="tomato", lw=1.4, label="Q_heat")
    style_axes(ax, title="Heating power", ylabel="kW")
    ax.legend(fontsize=9)

    ax = axes[1, 1]
    ax.plot(t_h, r["Q_cool"] / 1000, color="dodgerblue", lw=1.4, label="Q_cool")
    style_axes(ax, title="Cooling power", ylabel="kW")
    ax.legend(fontsize=9)

    ax = axes[2, 1]
    ax.plot(t_h, r["Q_vent"] / 1000, color="seagreen", lw=1.4, label="Q_vent")
    ax.axhline(0, color="gray", lw=0.8, ls="--")
    style_axes(ax, title="Dead band ventilation", ylabel="kW")
    ax.legend(fontsize=9)

    # Steady-state annotation
    UA_static  = UA_FACADE_W_K + UA_TUN_WALL_W_K + UA_SOIL_W_K
    T_tun_mean = r["T_ext"].mean() + T_TUN_OFFSET_DAY_C
    T_eq = (
        (UA_FACADE_W_K + UA_TUN_WALL_W_K) * T_tun_mean
        + UA_SOIL_W_K * T_SOIL_C
        + r["Q_int"].mean()
        + RHO_CP_AIR_J_M3_K * Q_STAIR_M3S * r["T_ext"].mean()
    ) / (UA_static + RHO_CP_AIR_J_M3_K * Q_STAIR_M3S)

    lines = [
        "Steady-state (no HVAC, mean):",
        f"  UA_fac+wall = {UA_FACADE_W_K + UA_TUN_WALL_W_K:.0f} W/K",
        f"  UA_soil     = {UA_SOIL_W_K:.0f} W/K",
        f"  T_tun_mean  = {T_tun_mean:.1f} °C",
        f"  T_soil      = {T_SOIL_C:.1f} °C",
        f"  Q_int_mean  = {r['Q_int'].mean()/1000:.1f} kW",
        f"  → T_eq ≈ {T_eq:.1f} °C",
        "",
        f"  C = {C_TOTAL_J_K:.2e} J/K",
        f"  τ = {C_TOTAL_J_K / UA_static / 3600:.0f} h",
    ]

    axes[3, 1].set_visible(False)
    ax_note = fig.add_axes([0.53, 0.02, 0.42, 0.22])
    ax_note.axis("off")
    ax_note.text(0.05, 0.95, "\n".join(lines), transform=ax_note.transAxes,
                 va="top", ha="left", fontsize=17, fontfamily="monospace",
                 bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Saved: {filename}")


# -----------------------------------------------------------------------------
# 3. SUMMARY PRINTOUT
# -----------------------------------------------------------------------------

def print_summary(r: dict):
    """Print annual energy, CO₂, cost, and monthly breakdown."""
    em = r["em"]
    print(f"T_in:  {r['T_in'].min():.1f}–{r['T_in'].max():.1f} °C")
    print(f"T_ext: {r['T_ext'].min():.1f}–{r['T_ext'].max():.1f} °C")
    print(f"Q_int: {r['Q_int'].min()/1e3:.1f}–{r['Q_int'].max()/1e3:.1f} kW")

    print(f"\n=== Annual Energy & Emissions ===")
    print(f"E_heating : {em['E_heat_total_kWh']:>8.0f} kWh")
    print(f"E_cooling : {em['E_cool_total_kWh']:>8.0f} kWh")
    print(f"E_fans    : {em['E_fan_total_kWh']:>8.0f} kWh")
    print(f"E_curtain : {em['E_curtain_total_kWh']:>8.0f} kWh")
    print(f"E_TOTAL   : {em['E_annual_kWh']:>8.0f} kWh")
    print(f"CO2_TOTAL : {em['CO2_annual_kgCO2']:>8.0f} kgCO₂")
    print(f"Curtain h : {r['curtain_on'].sum():>8} h")

    monthly = pd.DataFrame({
        "E_heat_kWh":    em["E_heat_kWh"],
        "E_cool_kWh":    em["E_cool_kWh"],
        "E_fan_kWh":     em["E_fan_kWh"],
        "E_curtain_kWh": em["E_curtain_kWh"],
        "CO2_total_kg":  em["CO2_total_kg"],
        "intensity":     r["co2_intensity"],
    }, index=r["dates"])

    ms = monthly.resample("ME").sum()
    ms["intensity_mean"] = monthly.resample("ME")["intensity"].mean()

    print("\n=== Monthly breakdown ===")
    print(ms[["E_heat_kWh", "E_cool_kWh", "E_fan_kWh", "E_curtain_kWh", "CO2_total_kg", "intensity_mean"]].round(0).to_string())

    print(f"\n=== Annual Cost ===")
    print(f"Heating : {em['cost_heat_eur']:>8.0f} €")
    print(f"Cooling : {em['cost_cool_eur']:>8.0f} €")
    print(f"Fans    : {em['cost_fan_eur']:>8.0f} €")
    print(f"Curtain : {em['cost_curtain_eur']:>8.0f} €")
    print(f"TOTAL   : {em['cost_annual_eur']:>8.0f} €")

    # Comfort hours
    above_27 = (r["T_in"] > 27.0).sum()
    above_26 = (r["T_in"] > 26.0).sum()
    below_18 = (r["T_in"] < 18.0).sum()
    total = len(r["T_in"])
    print(f"\n=== Comfort ===")
    print(f"Hours T_in > 27°C: {above_27:>5} / {total}  ({above_27/total*100:.1f}%)")
    print(f"Hours T_in > 26°C: {above_26:>5} / {total}  ({above_26/total*100:.1f}%)")
    print(f"Hours T_in < 18°C: {below_18:>5} / {total}  ({below_18/total*100:.1f}%)")

    is_service = np.array([not (1 <= ts.hour < 5) for ts in r["dates"]])

    above_26_service = ((r["T_in"] > 26.0) & is_service).sum()
    below_18_service = ((r["T_in"] < 18.0) & is_service).sum()
    below_18_night   = ((r["T_in"] < 18.0) & ~is_service).sum()
    service_hours    = is_service.sum()

    print(f"Service hours      : {service_hours} / {total}")
    print(f"Hours T_in > 26°C (service): {above_26_service} ({above_26_service/service_hours*100:.1f}%)")
    print(f"Hours T_in < 18°C (service): {below_18_service} ({below_18_service/service_hours*100:.1f}%)")
    print(f"Hours T_in < 18°C (night)  : {below_18_night} ({below_18_night/4*8784/8784:.0f}h expected max {8784//5})")


# -----------------------------------------------------------------------------
# 4. MAIN
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    r = run_simulation()
    plot_results(r)
    print_summary(r)
    import collections
    is_service = np.array([not (1 <= ts.hour < 5) for ts in r["dates"]])
    cold_service = [r["dates"][i].hour for i in range(len(r["T_in"]))
                    if r["T_in"][i] < 18.0 and is_service[i]]
    print(collections.Counter(sorted(cold_service)))