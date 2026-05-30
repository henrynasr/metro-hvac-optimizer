# =============================================================================
# simulation.py — Full simulation pipeline, one platform zone
# run_simulation()    → ODE + HVAC + emissions + humidity + latent + comfort
# classify_comfort()  → three-tier T / RH / combined
# print_summary()     → structured console output
# plot_results()      → 4×2 panel plot (to be redone)
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
    Q_STAIR_M3S, COP_COOL, ELEC_PRICE_EUR_KWH,
)
from occupancy import load_profiles, build_Q_array
from regulation import dT_dt, build_Q_hvac_array
from utils import load_data, style_axes
from emissions import load_co2_intensity, compute_emissions
from humidity import compute_humidity


# -----------------------------------------------------------------------------
# 1. COMFORT THRESHOLDS
# -----------------------------------------------------------------------------

T_COMFORT_LOW   = 18.0    # °C
T_COMFORT_HIGH  = 26.0
T_MILD_LOW      = 14.0
T_MILD_HIGH     = 28.0

RH_COMFORT_LOW  = 0.40
RH_COMFORT_HIGH = 0.60
RH_MILD_LOW     = 0.35
RH_MILD_HIGH    = 0.65


# -----------------------------------------------------------------------------
# 2. COMFORT CLASSIFICATION
# -----------------------------------------------------------------------------

def classify_comfort(
    T_in: np.ndarray,
    RH_in: np.ndarray,
    dates: pd.DatetimeIndex,
) -> dict:
    """
    Three-tier comfort classification on service hours (excludes 01h–05h).

    Temperature:
        comfort     18–26°C
        mild        14–18°C  or  26–28°C
        discomfort  <14°C    or  >28°C

    Humidity:
        comfort     40–60%
        mild        35–40%   or  60–65%
        discomfort  <35%     or  >65%

    Combined:
        comfort     T comfort  AND  RH comfort
        discomfort  T discomfort  OR  RH discomfort
        mild        everything else
    """
    is_service = np.array([not (1 <= ts.hour < 5) for ts in dates])
    n_service = is_service.sum()

    # --- Temperature ---
    t_comfort    = (T_in >= T_COMFORT_LOW) & (T_in <= T_COMFORT_HIGH)
    t_mild       = (((T_in >= T_MILD_LOW) & (T_in < T_COMFORT_LOW))
                  | ((T_in > T_COMFORT_HIGH) & (T_in <= T_MILD_HIGH)))
    t_discomfort = (T_in < T_MILD_LOW) | (T_in > T_MILD_HIGH)

    # --- Humidity ---
    rh_comfort    = (RH_in >= RH_COMFORT_LOW) & (RH_in <= RH_COMFORT_HIGH)
    rh_mild       = (((RH_in >= RH_MILD_LOW) & (RH_in < RH_COMFORT_LOW))
                   | ((RH_in > RH_COMFORT_HIGH) & (RH_in <= RH_MILD_HIGH)))
    rh_discomfort = (RH_in < RH_MILD_LOW) | (RH_in > RH_MILD_HIGH)

    # --- Combined ---
    combined_comfort    = t_comfort & rh_comfort
    combined_discomfort = t_discomfort | rh_discomfort
    combined_mild       = ~combined_comfort & ~combined_discomfort

    def pct(mask):
        return (mask & is_service).sum() / n_service * 100.0

    return {
        # Service hours info
        "is_service": is_service,
        "n_service": n_service,

        # Temperature (% of service hours)
        "T_comfort_pct":    pct(t_comfort),
        "T_mild_pct":       pct(t_mild),
        "T_discomfort_pct": pct(t_discomfort),

        # Humidity (% of service hours)
        "RH_comfort_pct":    pct(rh_comfort),
        "RH_mild_pct":       pct(rh_mild),
        "RH_discomfort_pct": pct(rh_discomfort),

        # Combined (% of service hours)
        "combined_comfort_pct":    pct(combined_comfort),
        "combined_mild_pct":       pct(combined_mild),
        "combined_discomfort_pct": pct(combined_discomfort),

        # Boolean arrays (for downstream analysis / monthly breakdowns)
        "T_comfort": t_comfort, "T_mild": t_mild, "T_discomfort": t_discomfort,
        "RH_comfort": rh_comfort, "RH_mild": rh_mild, "RH_discomfort": rh_discomfort,
        "combined_comfort": combined_comfort, "combined_mild": combined_mild,
        "combined_discomfort": combined_discomfort,
    }


# -----------------------------------------------------------------------------
# 3. SIMULATION ENGINE
# -----------------------------------------------------------------------------

def run_simulation(
    sim_start: str = "2024-01-01",
    sim_end: str = "2024-12-31",
    T0: float = 15.0,
    weather_csv: str = "data/raw/paris_weather.csv",
    xlsx_path: str = "data/raw/Defense_Occupation_Normalised.xlsx",
    rte_path: str = "data/raw/eco2mix_2024.csv",
) -> dict:
    """
    Full pipeline: ODE → HVAC decomposition → emissions → humidity → latent
    correction → comfort classification. Returns dict with all arrays and
    summary scalars. Every caller (Pareto, Sobol, standalone) gets the same
    complete result.
    """
    # --- Weather + occupancy inputs ---
    df       = load_data(weather_csv)
    df_sim   = df.loc[sim_start:sim_end].copy()
    dates    = df_sim.index
    T_ext    = df_sim["temperature_2m"].values
    RH_ext   = df_sim["relative_humidity_2m"].values / 100.0
    t_array  = np.arange(len(dates)) * 3600.0

    profiles        = load_profiles(xlsx_path)
    Q_int, n_people = build_Q_array(dates, profiles)

    # --- ODE solve ---
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

    # --- T_tun per hour ---
    T_tun = np.array([
        T_ext[i] + (T_TUN_OFFSET_NIGHT_C if 1 <= dates[i].hour < 5 else T_TUN_OFFSET_DAY_C)
        for i in range(len(dates))
    ])

    # --- Post-hoc HVAC decomposition ---
    (Q_total, Q_heat, Q_cool, Q_vent,
     Q_water_heat, Q_water_cool, T_hw, T_cw,
     Q_air_m3s_arr, curtain_on) = build_Q_hvac_array(T_in, T_ext, n_people, Q_int, dates)

    # --- Emissions ---
    co2_intensity = load_co2_intensity(rte_path, dates)
    em = compute_emissions(Q_heat, Q_cool, Q_vent, Q_air_m3s_arr, curtain_on,
                           co2_intensity, T_ext, T_hw)

    # --- Build results dict (needed by compute_humidity) ---
    r = {
        "dates": dates, "t_array": t_array,
        "T_ext": T_ext, "T_in": T_in, "T_tun": T_tun,
        "RH_ext": RH_ext,
        "Q_int": Q_int, "n_people": n_people,
        "Q_total": Q_total, "Q_heat": Q_heat, "Q_cool": Q_cool, "Q_vent": Q_vent,
        "Q_water_heat": Q_water_heat, "Q_water_cool": Q_water_cool,
        "T_hw": T_hw, "T_cw": T_cw,
        "Q_air_m3s_arr": Q_air_m3s_arr,
        "curtain_on": curtain_on,
        "co2_intensity": co2_intensity,
        "em": em,
    }

    # --- Humidity layer ---
    humidity = compute_humidity(r)
    r["humidity"] = humidity

    # --- Latent cooling correction ---
    P_latent_elec = humidity["latent_cool_W"] / COP_COOL          # W electrical per timestep
    E_latent_hourly = P_latent_elec / 1000.0                      # kWh (dt = 1h)
    E_latent_total = E_latent_hourly.sum()
    CO2_latent = (E_latent_hourly * co2_intensity / 1000.0).sum() # kgCO₂

    em["E_cool_latent_kWh"]  = E_latent_total
    em["E_cool_total_kWh"]  += E_latent_total
    em["E_annual_kWh"]      += E_latent_total
    em["CO2_annual_kgCO2"]  += CO2_latent
    em["cost_cool_eur"]     += E_latent_total * ELEC_PRICE_EUR_KWH
    em["cost_annual_eur"]   += E_latent_total * ELEC_PRICE_EUR_KWH

    # --- Comfort classification ---
    r["comfort"] = classify_comfort(T_in, humidity["RH_in"], dates)

    return r


# -----------------------------------------------------------------------------
# 4. SUMMARY PRINTOUT
# -----------------------------------------------------------------------------

def print_summary(r: dict):
    """Structured console output: energy, cost, comfort tiers, humidity."""
    em = r["em"]
    c  = r["comfort"]
    h  = r["humidity"]

    # --- Ranges ---
    print(f"T_in:  {r['T_in'].min():.1f} – {r['T_in'].max():.1f} °C")
    print(f"T_ext: {r['T_ext'].min():.1f} – {r['T_ext'].max():.1f} °C")
    print(f"Q_int: {r['Q_int'].min()/1e3:.1f} – {r['Q_int'].max()/1e3:.1f} kW")

    # --- Annual energy ---
    print(f"\n{'='*40}")
    print(f" Annual Energy & Emissions")
    print(f"{'='*40}")
    print(f"  Heating     {em['E_heat_total_kWh']:>8.0f} kWh")
    print(f"  Cooling     {em['E_cool_total_kWh']:>8.0f} kWh")
    print(f"  Fans        {em['E_fan_total_kWh']:>8.0f} kWh")
    print(f"  Curtain     {em['E_curtain_total_kWh']:>8.0f} kWh")
    print(f"  TOTAL       {em['E_annual_kWh']:>8.0f} kWh")
    print(f"  CO₂         {em['CO2_annual_kgCO2']:>8.0f} kgCO₂")

    # --- Annual cost ---
    print(f"\n{'='*40}")
    print(f" Annual Cost")
    print(f"{'='*40}")
    print(f"  Heating     {em['cost_heat_eur']:>8.0f} €")
    print(f"  Cooling     {em['cost_cool_eur']:>8.0f} €")
    print(f"  Fans        {em['cost_fan_eur']:>8.0f} €")
    print(f"  Curtain     {em['cost_curtain_eur']:>8.0f} €")
    print(f"  TOTAL       {em['cost_annual_eur']:>8.0f} €")

    # --- Comfort: Temperature ---
    print(f"\n{'='*40}")
    print(f" Comfort — Temperature (service hours)")
    print(f"{'='*40}")
    print(f"  Comfort   (18–26°C)    {c['T_comfort_pct']:>5.1f}%")
    print(f"  Mild      (14–18/26–28°C) {c['T_mild_pct']:>5.1f}%")
    print(f"  Discomfort (<14/>28°C) {c['T_discomfort_pct']:>5.1f}%")

    # --- Comfort: Humidity ---
    print(f"\n{'='*40}")
    print(f" Comfort — Humidity (service hours)")
    print(f"{'='*40}")
    print(f"  Comfort   (40–60%)     {c['RH_comfort_pct']:>5.1f}%")
    print(f"  Mild      (35–40/60–65%)  {c['RH_mild_pct']:>5.1f}%")
    print(f"  Discomfort (<35/>65%)  {c['RH_discomfort_pct']:>5.1f}%")

    # --- Comfort: Combined ---
    print(f"\n{'='*40}")
    print(f" Comfort — Combined (service hours)")
    print(f"{'='*40}")
    print(f"  Comfort   (both OK)    {c['combined_comfort_pct']:>5.1f}%")
    print(f"  Mild      (tolerable)  {c['combined_mild_pct']:>5.1f}%")
    print(f"  Discomfort (≥1 bad)    {c['combined_discomfort_pct']:>5.1f}%")

    # --- Humidity details ---
    print(f"\n{'='*40}")
    print(f" Humidity Details")
    print(f"{'='*40}")
    print(f"  Condensation hours     {h['hours_condensation']:>5}")
    print(f"  Latent cooling load    {h['latent_cool_total_kWh']:>8.0f} kWh thermal")
    print(f"  Latent elec (÷COP)    {em.get('E_cool_latent_kWh', 0):>8.0f} kWh elec")
    print(f"  Curtain runtime        {r['curtain_on'].sum():>5} h")

    # --- Monthly breakdown ---
    monthly = pd.DataFrame({
        "E_heat":    em["E_heat_kWh"],
        "E_cool":    em["E_cool_kWh"],
        "E_fan":     em["E_fan_kWh"],
        "E_curtain": em["E_curtain_kWh"],
        "CO2":       em["CO2_total_kg"],
    }, index=r["dates"])

    ms = monthly.resample("ME").sum()
    print(f"\n{'='*40}")
    print(f" Monthly Breakdown (kWh / kgCO₂)")
    print(f"{'='*40}")
    print(ms.round(0).to_string())


# -----------------------------------------------------------------------------
# 5. PLOT (unchanged — to be redone)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 5a. TEMPERATURES
# -----------------------------------------------------------------------------
 
def plot_temperatures(r: dict, filename: str = "images/temperatures.png"):
    """T_ext, T_in, T_tun, T_soil with comfort bands."""
    fig, ax = plt.subplots(figsize=(16, 6))
    dates = r["dates"]
 
    # Comfort bands
    ax.axhspan(T_COMFORT_LOW, T_COMFORT_HIGH, color="#27ae60", alpha=0.08)
    ax.axhspan(T_MILD_LOW, T_COMFORT_LOW, color="#f39c12", alpha=0.06)
    ax.axhspan(T_COMFORT_HIGH, T_MILD_HIGH, color="#f39c12", alpha=0.06)
 
    ax.plot(dates, r["T_ext"], color="steelblue", lw=0.4, alpha=0.5, label="T_ext")
    ax.plot(dates, r["T_tun"], color="darkorange", lw=0.3, alpha=0.3, label="T_tun")
    ax.axhline(T_SOIL_C, color="sienna", lw=1.0, ls=":", alpha=0.7,
               label=f"T_soil = {T_SOIL_C}°C")
    ax.plot(dates, r["T_in"], color="crimson", lw=0.6, label="T_in")
 
    # Band annotations (right edge)
    x_end = dates[-1]
    ax.annotate("Comfort\n18–26°C", xy=(x_end, 22), fontsize=7,
                color="#27ae60", fontweight="bold", ha="left",
                xytext=(10, 0), textcoords="offset points")
    ax.annotate("Mild", xy=(x_end, 16), fontsize=7,
                color="#b8860b", ha="left",
                xytext=(10, 0), textcoords="offset points")
    ax.annotate("Mild", xy=(x_end, 27), fontsize=7,
                color="#b8860b", ha="left",
                xytext=(10, 0), textcoords="offset points")
 
    style_axes(ax, title="Temperatures — Full Year 2024", ylabel="°C")
    ax.legend(loc="upper left", fontsize=9, framealpha=0.9)
 
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Saved: {filename}")
 
 
# -----------------------------------------------------------------------------
# 5b. HVAC DECOMPOSITION
# -----------------------------------------------------------------------------
 
def plot_hvac(r: dict, filename: str = "images/hvac.png"):
    """2×2 panels: heating, cooling, ventilation, headcount + Q_int."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=True)
    dates = r["dates"]
 
    ax = axes[0, 0]
    ax.fill_between(dates, r["Q_heat"] / 1000, color="tomato", alpha=0.6)
    ax.plot(dates, r["Q_heat"] / 1000, color="tomato", lw=0.3)
    style_axes(ax, title="Heating power", ylabel="kW")
 
    ax = axes[0, 1]
    ax.fill_between(dates, r["Q_cool"] / 1000, color="dodgerblue", alpha=0.6)
    ax.plot(dates, r["Q_cool"] / 1000, color="dodgerblue", lw=0.3)
    style_axes(ax, title="Cooling power", ylabel="kW")
 
    ax = axes[1, 0]
    ax.fill_between(dates, r["Q_vent"] / 1000, color="seagreen", alpha=0.5)
    ax.plot(dates, r["Q_vent"] / 1000, color="seagreen", lw=0.3)
    ax.axhline(0, color="gray", lw=0.6, ls="--")
    style_axes(ax, title="Dead band ventilation", ylabel="kW")
 
    ax = axes[1, 1]
    ax.fill_between(dates, r["n_people"], color="teal", alpha=0.4)
    ax.plot(dates, r["n_people"], color="teal", lw=0.3)
    ax2 = ax.twinx()
    ax2.plot(dates, r["Q_int"] / 1000, color="darkorchid", lw=0.5, alpha=0.7)
    ax2.set_ylabel("Q_int [kW]", fontsize=10, color="darkorchid")
    ax2.tick_params(axis="y", labelcolor="darkorchid")
    style_axes(ax, title="Headcount + internal gains", ylabel="Persons")
 
    fig.suptitle("HVAC Decomposition — Full Year 2024", fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {filename}")
 
 
# -----------------------------------------------------------------------------
# 5c. HUMIDITY
# -----------------------------------------------------------------------------
 
def plot_humidity(r: dict, filename: str = "images/humidity.png"):
    """RH_in over time with comfort bands."""
    fig, ax = plt.subplots(figsize=(16, 5))
    dates = r["dates"]
    RH_in = r["humidity"]["RH_in"] * 100  # 0–1 → %
 
    # Comfort bands
    ax.axhspan(RH_COMFORT_LOW * 100, RH_COMFORT_HIGH * 100,
               color="#27ae60", alpha=0.08)
    ax.axhspan(RH_MILD_LOW * 100, RH_COMFORT_LOW * 100,
               color="#f39c12", alpha=0.06)
    ax.axhspan(RH_COMFORT_HIGH * 100, RH_MILD_HIGH * 100,
               color="#f39c12", alpha=0.06)
 
    ax.plot(dates, RH_in, color="royalblue", lw=0.4, alpha=0.7)
 
    # Condensation markers
    cond = r["humidity"]["condensation"]
    if cond.any():
        ax.scatter(dates[cond], RH_in[cond], s=3, color="red", alpha=0.5,
                   zorder=5, label=f"Condensation ({cond.sum()} h)")
        ax.legend(loc="upper right", fontsize=9)
 
    style_axes(ax, title="Indoor Relative Humidity — Full Year 2024", ylabel="RH [%]")
    ax.set_ylim(10, 100)
 
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Saved: {filename}")
 
 
# -----------------------------------------------------------------------------
# 5d. WATER TEMPERATURES
# -----------------------------------------------------------------------------
 
def plot_water(r: dict, filename: str = "images/water_temps.png"):
    """Hot and cold water supply temps + T_ext context."""
    fig, ax = plt.subplots(figsize=(16, 5))
    dates = r["dates"]
 
    ax.plot(dates, r["T_ext"], color="steelblue", lw=0.3, alpha=0.4, label="T_ext")
    ax.plot(dates, r["T_hw"], color="tomato", lw=1.0, label="T_hw supply")
    ax.plot(dates, r["T_cw"], color="dodgerblue", lw=1.0, label="T_cw supply")
 
    # Reference lines
    ax.axhline(0, color="gray", lw=0.5, ls="--", alpha=0.5)
 
    style_axes(ax, title="Water Circuit Temperatures — Full Year 2024", ylabel="°C")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.9)
 
    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    plt.close()
    print(f"Saved: {filename}")
 
 
# -----------------------------------------------------------------------------
# 5e. ENERGY + CO₂ MONTHLY
# -----------------------------------------------------------------------------
 
def plot_energy_monthly(r: dict, filename: str = "images/energy_monthly.png"):
    """Stacked bar: monthly energy by component + CO₂ bars below."""
    em = r["em"]
 
    monthly = pd.DataFrame({
        "Heating":  em["E_heat_kWh"],
        "Cooling":  em["E_cool_kWh"],
        "Fans":     em["E_fan_kWh"],
        "Curtain":  em["E_curtain_kWh"],
        "CO2":      em["CO2_total_kg"],
    }, index=r["dates"])
    ms = monthly.resample("ME").sum()
 
    labels = [d.strftime("%b") for d in ms.index]
    x = np.arange(len(labels))
    w = 0.6
 
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True,
                                    gridspec_kw={"height_ratios": [3, 1]})
 
    # Energy stacked bar
    colors = {"Heating": "#e74c3c", "Cooling": "#3498db",
              "Fans": "#2ecc71", "Curtain": "#e67e22"}
    bottom = np.zeros(len(x))
    for col in ["Heating", "Cooling", "Fans", "Curtain"]:
        vals = ms[col].values
        ax1.bar(x, vals, w, bottom=bottom, label=col, color=colors[col], alpha=0.85)
        bottom += vals
 
    style_axes(ax1, title="Monthly Energy Breakdown", ylabel="kWh")
    ax1.legend(loc="upper left", fontsize=9, framealpha=0.9)
 
    # CO₂ bars
    ax2.bar(x, ms["CO2"].values, w, color="dimgray", alpha=0.8)
    style_axes(ax2, ylabel="kgCO₂")
    ax2.set_xticks(x)
    ax2.set_xticklabels(labels, fontsize=10)
 
    fig.suptitle("Energy & CO₂ — Monthly 2024", fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {filename}")
 
 
# -----------------------------------------------------------------------------
# 5f. COMFORT TIERS MONTHLY
# -----------------------------------------------------------------------------
 
def _monthly_tier_pct(
    tier_arr: np.ndarray,
    is_service: np.ndarray,
    dates: pd.DatetimeIndex,
) -> np.ndarray:
    """Percentage of service hours in a given tier, per month."""
    df = pd.DataFrame({
        "tier": tier_arr & is_service,
        "service": is_service,
    }, index=dates)
    m_tier = df["tier"].resample("ME").sum()
    m_serv = df["service"].resample("ME").sum()
    return (m_tier / m_serv * 100).values
 
 
def plot_comfort_monthly(r: dict, filename: str = "images/comfort_monthly.png"):
    """3 stacked bar panels: temperature, humidity, combined comfort tiers."""
    c = r["comfort"]
    dates = r["dates"]
    is_service = c["is_service"]
 
    labels = [d.strftime("%b") for d in
              pd.date_range(dates[0], dates[-1], freq="ME")]
    x = np.arange(len(labels))
    w = 0.6
 
    colors = {"comfort": "#27ae60", "mild": "#f1c40f", "discomfort": "#e74c3c"}
 
    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
 
    # --- Temperature ---
    t_c = _monthly_tier_pct(c["T_comfort"], is_service, dates)
    t_m = _monthly_tier_pct(c["T_mild"], is_service, dates)
    t_d = _monthly_tier_pct(c["T_discomfort"], is_service, dates)
 
    axes[0].bar(x, t_c, w, label="Comfort (18–26°C)", color=colors["comfort"], alpha=0.85)
    axes[0].bar(x, t_m, w, bottom=t_c, label="Mild (14–18 / 26–28°C)", color=colors["mild"], alpha=0.85)
    axes[0].bar(x, t_d, w, bottom=t_c + t_m, label="Discomfort (<14 / >28°C)", color=colors["discomfort"], alpha=0.85)
    style_axes(axes[0], title="Temperature Comfort", ylabel="%")
    axes[0].set_ylim(0, 105)
    axes[0].legend(loc="lower right", fontsize=8, framealpha=0.9)
 
    # --- Humidity ---
    rh_c = _monthly_tier_pct(c["RH_comfort"], is_service, dates)
    rh_m = _monthly_tier_pct(c["RH_mild"], is_service, dates)
    rh_d = _monthly_tier_pct(c["RH_discomfort"], is_service, dates)
 
    axes[1].bar(x, rh_c, w, label="Comfort (40–60%)", color=colors["comfort"], alpha=0.85)
    axes[1].bar(x, rh_m, w, bottom=rh_c, label="Mild (35–40 / 60–65%)", color=colors["mild"], alpha=0.85)
    axes[1].bar(x, rh_d, w, bottom=rh_c + rh_m, label="Discomfort (<35 / >65%)", color=colors["discomfort"], alpha=0.85)
    style_axes(axes[1], title="Humidity Comfort", ylabel="%")
    axes[1].set_ylim(0, 105)
    axes[1].legend(loc="lower right", fontsize=8, framealpha=0.9)
 
    # --- Combined ---
    cb_c = _monthly_tier_pct(c["combined_comfort"], is_service, dates)
    cb_m = _monthly_tier_pct(c["combined_mild"], is_service, dates)
    cb_d = _monthly_tier_pct(c["combined_discomfort"], is_service, dates)
 
    axes[2].bar(x, cb_c, w, label="Comfort (both OK)", color=colors["comfort"], alpha=0.85)
    axes[2].bar(x, cb_m, w, bottom=cb_c, label="Mild (tolerable)", color=colors["mild"], alpha=0.85)
    axes[2].bar(x, cb_d, w, bottom=cb_c + cb_m, label="Discomfort (≥1 bad)", color=colors["discomfort"], alpha=0.85)
    style_axes(axes[2], title="Combined Comfort", ylabel="%")
    axes[2].set_ylim(0, 105)
    axes[2].legend(loc="lower right", fontsize=8, framealpha=0.9)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels(labels, fontsize=10)
 
    fig.suptitle("Comfort Tiers — Monthly 2024 (service hours only)", fontsize=14, y=1.01)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved: {filename}")
 
 
# -----------------------------------------------------------------------------
# 5g. DISPATCHER
# -----------------------------------------------------------------------------
 
def plot_results(r: dict):
    """Generate all plots."""
    plot_temperatures(r)
    plot_hvac(r)
    plot_humidity(r)
    plot_water(r)
    plot_energy_monthly(r)
    plot_comfort_monthly(r)


# -----------------------------------------------------------------------------
# 6. MAIN
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    r = run_simulation()
    plot_results(r)
    print_summary(r)