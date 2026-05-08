# Energy Twin

HVAC sensitivity & optimization toolkit for metro stations. Public data only.

## What this is

A Python toolkit for HVAC sensitivity analysis and regulation strategy optimization on metro stations, built exclusively on public data (Open-Meteo, IDFM/RATP, RTE, ADEME). Longer-term goal: a transfer learning methodology — calibrate a thermal + regulation model on a data-rich reference station, then deploy it on a new station with zero operational history.

Current state: 30-year Paris climate analysis, a physically-grounded lumped-capacitance thermal model with a full regulation layer, sourced and documented parameters, and a quantitative Sobol sensitivity analysis on envelope and load parameters.

---

## Model architecture

One platform zone (one side, one AHU). Two symmetric zones per station — identical, not thermally coupled. Model one, apply to both.

**ODE (Session 11):**

```
C · dT_in/dt = (UA_f + ρcp·V̇_inf)·(T_tun − T_in)
             +  UA_s·(T_soil − T_in)
             +  Q_int
             −  Q_hvac
```

T_ext is not a direct envelope boundary — the station is fully buried. Two boundaries:
- **Facade (PSD glass wall)** → tunnel air: T_tun = T_ext + 5°C
- **Structural box (concrete + soil)** → stable ground: T_soil = 15°C

Infiltration V̇_inf is dynamic: computed from train frequency by hour and day-type (peak 1.155 m³/s, off-peak 0.578 m³/s, night 0).

**Key parameters (all sourced — see docs/parameters.md):**

| Parameter | Value | Source |
|---|---|---|
| Platform zone | 55m × 4m × 4.2m | MP14 rolling stock (DUP 2012) + Fruin LOS C |
| U_facade | 7.0 W/m²K | ASHRAE 90.1 + EN ISO 10077 |
| U_soil | 0.41 W/m²K | ISO 13370 + ASHRAE Ch.25 |
| C_total | 1.73×10⁸ J/K | Air + concrete surface layer (12 cm effective) |
| Peak occupancy | 250 persons/side | Fruin LOS C, Metrolinx DS-12 (2024) |
| Sensible heat | 75 W/person | ASHRAE Fundamentals 2013, Ch.18, Table 1 |
| Baseline equipment | 5 kW | SEAM4US Barcelona × LED factor |

---

## Findings

### Paris climate — 30 years (1996–2025)

- **Warming trend: +0.48 °C/decade**, p = 0.0002. Consistent with Western European trends.
- Diurnal swing: ~6.7 °C average. Min at 5–6 AM, max at 2–3 PM.
- Temperature vs humidity: Pearson r = −0.58. Clausius-Clapeyron directly visible — relevant to HVAC dehumidification load.

### Sobol sensitivity analysis — 6 parameters, 2 metrics

Saltelli sampling, N=512, 7168 ODE solves. Inputs: UA, C, watts/person, baseline load, peak headcount, T0. Metrics: peak T_in, % hours > 26°C.

- **UA dominates both metrics** (ST = 0.43 / 0.59). ~2× any other input.
- **Model is additive** — no significant pairwise interactions.
- **C matters for peak T_in but not for hot hours** (ST = 0.15 vs 0.05).
- **T0 is noise** — initial conditions fade in a few time constants.

Next Sobol run will target the new parameters: T_tun offset, infiltration η, U_facade, U_soil — the dominant unknowns in the updated ODE.

### Regulation layer

Five-zone setpoint law driven by T_ext:

| T_ext | Mode | T_in target |
|---|---|---|
| < −1°C | Anti-freeze | ≥ 5°C |
| −1 to 6°C | Heating (proportional) | T_ext + 6°C |
| 6 to 12°C | Heating (fixed) | 12°C |
| 12 to 26°C | Dead band | — (ventilation only) |
| 26 to 31°C | Cooling (fixed) | 26°C |
| > 31°C | Cooling (proportional) | T_ext − 5°C |

Airflow: Q_total = Q_overpressure (2500 m³/h, constant) + Q_occupancy (25 m³/h/person). Range 2500–9625 m³/h. Overpressure and hygiene demands are independent — they add, not max.

HVAC power uses fixed T_blow: 15°C in cooling mode, 30°C in heating mode. Q_hvac = ρcp × Q_air × (T_in − T_blow), with shut-off when T_in already on the right side of T_set.

Power decomposed into Q_heat / Q_cool / Q_vent — mutually exclusive, readable at a glance.

**Observed July 2024:** T_in 21–25°C while T_ext 10–24°C. Dead-band ventilation (Q_vent 5–17 kW) does most of the work — outdoor air cooler than the zone pulls T_in down continuously. Active cooling fires only on the hottest days.

---

## Data

- **Open-Meteo Historical Weather API** — ERA5 reanalysis, 30 years hourly Paris (1996–2025). Raw CSV gitignored.
- **RATP — Fréquentation du pôle La Défense** (IDFM open data). Normalised hourly profiles by day-type (JOHV/JOVS/WKD), anchored to JOHV 18h = 250 persons. File: `data/raw/Defense_Occupation_Normalised.xlsx`.

---

## Scripts

- `fetch_weather.py` — pulls 30y Paris weather from Open-Meteo → `data/raw/paris_weather.csv`
- `inspect_weather.py` — shape, dtypes, summary stats, missing values
- `plot_weather.py` — 5 climate analysis plots → `images/`
- `constants.py` — single source of truth for all model parameters (sourced + flagged)
- `occupancy.py` — RATP profiles, day-type dispatch, `build_Q_array`, `v_inf_m3s`
- `regulation.py` — `T_setpoint`, `airflow_total`, `dT_dt`, `build_Q_hvac_array`
- `utils.py` — `load_data`, `style_axes`
- `thermal_model.py` — single-run simulation, 4×2 panel plot
- `sweep.py` — 20×20 UA × w/p heatmap (pending update to new ODE)
- `sobol.py` — 6-parameter Sobol via SALib (pending update to new ODE)
- `docs/parameters.md` — full parameter table with sources, derivations, Sobol flags

---

## Setup

```bash
git clone https://github.com/henrynasr/energy-twin
cd energy-twin
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
python fetch_weather.py         # run first — pulls weather data
python thermal_model.py
```

---

## Status

**Week 2, Session 11.** ODE upgraded to split boundary conditions (facade→T_tun, structure→T_soil) with dynamic infiltration V̇_inf(hour, day_type). All parameters re-sourced and documented in docs/parameters.md. T_blow fixed at 15°C/30°C (cooling/heating). Airflow corrected to 250 persons/side, additive overpressure logic. Validated on July 2024 — physically consistent. Next: update sweep.py and sobol.py to new ODE, then Sobol on the new parameter set (T_tun offset, η, U_facade, U_soil).
