# Energy Twin

HVAC sensitivity & optimization toolkit for metro stations. Public data only.

## What this is

A Python toolkit for HVAC sensitivity analysis and regulation strategy optimization on metro stations, built exclusively on public data sources (Open-Meteo, Météo-France, IDFM/RATP, RTE, ADEME). The longer-term goal is a transfer learning methodology: calibrate a thermal + regulation model on a data-rich reference station, then deploy it on a new station with zero operational history.

Current state: weather data pipeline, 30-year Paris climate analysis, a lumped-capacitance thermal model driven by real RATP occupancy profiles, a full regulation layer (setpoint law, airflow modulation, HVAC power decomposition), and a quantitative Sobol sensitivity analysis on envelope and load parameters.

---

## Findings

### Paris climate — 30 years (1996–2025)

30 years of hourly ERA5 reanalysis data pulled from Open-Meteo. Key results:

- **Warming trend: +0.48 °C/decade**, p = 0.0002. Statistically significant, consistent with reported Western European trends.
- Diurnal swing: minimum at 5–6 AM (~8.7 °C), maximum at 2–3 PM (~15.4 °C), average swing ~6.7 °C.
- Temperature vs relative humidity: Pearson r = −0.58. Cold air clusters at high humidity; warm air spans wide RH range — a direct consequence of Clausius-Clapeyron (air water-holding capacity scales ~7%/°C). Relevant to HVAC dehumidification load: summer cooling is also water removal.
- Notable outliers: 2003 European heatwave (July peak ~23.5 °C monthly mean), cold winter of 2010, cluster of warm summers in the 2020s.

### Thermal model

Lumped-capacitance ODE: `C·dT_in/dt = UA·(T_ext − T_in) + Q_internal − Q_hvac`, driven by real Paris weather and real RATP occupancy profiles from Pôle La Défense.

Occupancy data pre-processing: filtered to post-COVID years, aggregated to hourly profiles by day-type (JOHV / JOVS / WKD), normalised on JOHV-18h = 100%, calibrated to 400 persons/side as peak. Day-type dispatch is holiday-aware (Zone C 2024 school calendar).

Three physical behaviors consistently observed across all simulation windows:
- **Offset:** T_in sits above T_ext at peak — driven by internal gains (occupancy + lighting + equipment).
- **Lag:** T_in peaks several hours after T_ext, consistent with thermal time constant τ = C/UA ≈ 2.8 h.
- **Damping:** T_in is markedly smoother than both T_ext and occupancy spikes — thermal mass acts as a low-pass filter.

Parameters (UA, C) are order-of-magnitude estimates, not calibrated against a real station.

### Sensitivity sweep — UA × watts per person

400 thermal model runs sweeping UA (envelope conductance, ±40% around baseline) and watts per person (70–120 W, seated to active occupancy). Peak T_in ranges from ~28 °C to ~35 °C across the parameter space — a 7 °C spread under realistic conditions for a buried, unventilated station in summer. The gradient is diagonal, meaning both parameters contribute independently to the peak response. Their relative weights are quantified in the Sobol analysis below.

### Sobol sensitivity analysis — 6 parameters, 2 metrics

Variance-based global sensitivity analysis (SALib, Saltelli sampling, N=512, 7168 ODE solves). Six inputs: UA, C, watts per person, baseline equipment load, peak headcount, initial temperature. Two metrics: peak T_in and % hours T_in > 26 °C.

Key findings:
- **UA dominates both metrics** (ST = 0.43 on peak, ST = 0.59 on hot-hour fraction). Envelope conductance is the single biggest lever — roughly 2× any other input.
- **Model is additive.** S1 ≈ ST for every parameter, all pairwise interaction indices statistically indistinguishable from zero.
- **C matters for peaks but not for hot hours** (ST = 0.15 vs 0.05). Thermal mass smooths the peak; it does not change how often the system runs hot once the load profile is fixed.
- **T0 is noise.** Initial conditions fade within a few thermal time constants and are irrelevant for weekly statistics.
- **wpp and PEOPLE_PEAK split the occupancy load** roughly evenly (~ST 0.18 each on peak).

This closes the visual sweep by ranking inputs quantitatively. The next sensitivity analysis will target regulation parameters, not envelope/load.

### Regulation layer

A physically-grounded HVAC regulation model built from first principles, covering a 2-AHU underground metro station.

**Setpoint law — 5 zones:**

| T_ext range | T_in target | Mode |
|---|---|---|
| < −1 °C | ≥ 5 °C | Anti-freeze |
| −1 to 6 °C | T_ext + 6 °C | Heating, proportional to outdoor |
| 6 to 12 °C | 12 °C | Heating, fixed target |
| 12 to 26 °C | — | Dead band, ventilation only |
| 26 to 31 °C | 26 °C | Cooling, fixed target |
| > 31 °C | T_ext − 5 °C | Cooling, proportional to outdoor |

**Airflow modulation:** proportional to instantaneous occupancy, bounded between 10,000 m³/h (minimum overpressure at night) and 22,000 m³/h (full capacity at 800-person peak, 2 AHUs × 11,000 m³/h). Sized on 25 m³/h/person at peak + 10% safety factor. In the absence of real-time CO₂ sensor data, a backup pre-defined hourly airflow schedule (pending integration) provides fallback control.

**HVAC power:** `Q_hvac = ρ × Q_air × Cp × ΔT_controlled`, where ΔT is clipped to [−12, +5] K — heating limited by the anti-freeze worst case (T_ext = −7 °C, T_set = 5 °C), cooling limited by the T_ext > 31 °C zone (T_set = T_ext − 5, ΔT = 5 K by construction).

**Power decomposition:** Q_hvac is split into three non-overlapping components — Q_heat (active heating), Q_cool (active cooling), Q_vent (dead-band ventilation heat exchange). When one is non-zero, the other two are zero. This makes the controller's behavior readable at a glance across seasons.

**Observed behavior:**
- *Summer dead band (July 2024):* Q_heat = 0, Q_cool fires briefly when T_ext spikes above 26 °C, Q_vent dominates as outdoor air at ~20 °C pulls T_in down from internal gains.
- *Winter heating (January 2024):* Q_heat is the base signal, negative (adding heat). Q_cool spikes during rush hours when internal gains (50 kW peak from occupancy) push T_in above T_set despite cold outdoor conditions — ventilation of cold outside air is mandatory for hygiene and overpressure, and that air is only heated to 12 °C, so if T_in reaches 15 °C the system correctly reads it as needing cooling. This is physically expected, not a model artifact.
- *Heatwave (August 2003):* Q_cool dominates, Q_heat near zero, Q_vent active in dead-band hours.

Next step: water regime law (glycol temperature as a function of T_ext, with hysteresis) to set a physically-grounded cap on Q_hvac_max.

---

## Data

- **Open-Meteo Historical Weather API** (ERA5 reanalysis), 30 years of Paris hourly weather (1996–2025). Variables: 2 m temperature, 2 m relative humidity. Raw CSV not committed (`.gitignore`).
- **RATP — Fréquentation du pôle La Défense** (Île-de-France Mobilités open data). Filtered to post-COVID years, aggregated to hourly profiles by day-type, normalised on JOHV-18h = 100%. Processed file tracked in the repo: `data/raw/Defense_Occupation_Normalised.xlsx`.

---

## Scripts

- `fetch_weather.py` — pulls 30 years of hourly Paris weather from Open-Meteo, saves to `data/raw/paris_weather.csv`.
- `inspect_weather.py` — loads the CSV, prints shape, dtypes, summary stats, missing values.
- `plot_weather.py` — generates climate analysis plots into `images/`.
- `occupancy.py` — reads RATP normalised profiles, dispatches day-types (JOHV/JOVS/WKD) with 2024 Zone C holiday calendar, builds `(Q_array, n_people_array)` from a datetime index.
- `regulation.py` — setpoint law (`T_setpoint`), airflow modulation (`airflow_total`), ODE slope function (`dT_dt`), and post-hoc HVAC power decomposition (`build_Q_hvac_array`).
- `utils.py` — shared functions: `load_data`, `style_axes`.
- `thermal_model.py` — single-run simulation: lumped-capacitance ODE driven by weather + occupancy + regulation. 4×2 panel plot (temperatures, internal load, headcount, total Q_hvac / heating, cooling, ventilation breakdown).
- `sweep.py` — 20×20 sensitivity sweep on (UA, watts_per_person), heatmap of peak T_in.
- `sobol.py` — Sobol global sensitivity analysis on 6 parameters via SALib, two output metrics.

---

## Setup

```bash
git clone https://github.com/henrynasr/energy-twin
cd energy-twin
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
python fetch_weather.py         # run this first — pulls the weather data
python thermal_model.py
python sweep.py
python sobol.py
```

---

## Status

**Week 2, Session 9.** Regulation layer complete: setpoint law, airflow modulation, HVAC power decomposition (heat/cool/vent). Validated on three regimes — summer dead band, winter heating, 2003 heatwave. Next: water regime law (glycol temperature vs T_ext with hysteresis) to cap controller power physically, then Sobol on regulation parameters.
