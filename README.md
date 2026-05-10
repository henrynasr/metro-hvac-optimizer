# Energy Twin — HVAC Sensitivity & Optimization Toolkit

Physics-based thermal model of a Paris underground metro stop — one side platform, one AHU.
Built on public data only. No proprietary files.

---

## What it does

Models the thermal dynamics of one buried platform zone using a lumped-capacitance ODE driven by real Paris weather (Open-Meteo ERA5) and real RATP occupancy data. The regulation layer replicates a real AHU control strategy: two-zone setpoint law, airflow modulation, water circuit scheduling, and hysteresis logic. A Sobol global sensitivity analysis identifies the parameters that matter most.

**ODE:**

```
C · dT_in/dt = (UA_facade + UA_tun_wall + ρcp·V̇_inf)·(T_tun − T_in)
             +  UA_soil·(T_soil − T_in)
             +  Q_int
             +  ρcp·Q_stair·(T_ext − T_in)
             −  Q_hvac
```

- `T_tun = T_ext + 10°C` during service (05h–23h), `T_ext + 5°C` at night — train braking, motor losses, passenger heat. [SOBOL 5–15 / 3–8]
- `T_soil = 15°C` — stable ground temperature (BRGM). Applied to closed concrete side only.
- `UA_tun_wall` — concrete wall above PSD glass (1.4m × 55m), faces tunnel. U = 2.5 W/m²K.
- `V̇_inf` — dynamic PSD infiltration, f(hour, day-type), exchange efficiency method.
- `Q_stair` — passive outdoor air via open stair entrance. Q = V_air × A_stair = 1.25 m³/s fixed.
- `Q_int` — sensible heat from occupancy (75 W/person, ASHRAE) + equipment (500 W LED/screens, corridor only).
- `Q_hvac` — AHU power. T_blow fixed: 15°C cooling, 30°C heating.

**Setpoint law:**

- T_ext < 15°C → heat to 21°C
- 15°C ≤ T_ext ≤ 20°C → dead band (ventilation only)
- T_ext > 20°C → cool to T_ext − 6°C (e.g. 34°C when T_ext = 40°C)

**Water regime:**

Hot circuit: 50→35°C supply over T_ext −7→12°C. Off above 12°C (restarts at 10°C).
Cold circuit: 12→8°C supply over T_ext 26→31°C. Off below 26°C (restarts at 28°C).
`Q_water = Q_air × ρcp_air × dT_air / (ρ_glycol × Cp_glycol × ΔT_water)` where `T_mix = 0.7·T_in + 0.3·T_ext`.

---

## Data

- **Open-Meteo Historical Weather API** — ERA5 reanalysis, 30 years hourly Paris (1996–2025). Raw CSV gitignored.
- **RATP — Fréquentation du pôle La Défense** (IDFM open data). Normalised hourly profiles by day-type (JOHV/JOVS/WKD), anchored to JOHV 18h = 250 persons.

---

## Scripts

| File | Role |
|---|---|
| `fetch_weather.py` | Pull 30y Paris weather from Open-Meteo |
| `constants.py` | Single source of truth — all parameters with units and Sobol ranges |
| `occupancy.py` | RATP profiles, day-type dispatch, `v_inf_m3s` |
| `regulation.py` | `T_setpoint`, `airflow_total`, `dT_dt`, `build_Q_hvac_array`, water regime |
| `thermal_model.py` | Single-run simulation, 4×2 panel plot |
| `sobol_A.py` | Sobol A (27 params, screen) + C (5 survivors, N=512) via SALib. Self-contained parametric ODE. |
| `sobol_B.py` | Sobol B — water regime sensitivity, post-hoc (no ODE per row). August week. |
| `utils.py` | `load_data`, `style_axes` |
| `docs/parameters.md` | Full parameter table — values, sources, derivations, Sobol ranges |

---

## Setup

```bash
git clone https://github.com/henrynasr/energy-twin
cd energy-twin
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
python fetch_weather.py         # run once — pulls weather data
python thermal_model.py
```

---

## Key findings

**August 2024:** T_in 24–32°C, T_ext 18–37°C. Active cooling fires on the 2 hottest days (peak ~40 kW). Dead-band ventilation active throughout. Steady-state T_eq ≈ 28.3°C without HVAC — confirms cooling is necessary in summer.

**January 2024:** T_in 10–21°C, T_ext −4–4°C. Heating fires at startup then soil + Q_int maintain T_in passively above setpoint for most of the week. Staircase infiltration pulls T_in down during cold spells.

**Sobol C (5 params, N=512, July week):**
- `T_TUN_OFFSET_DAY` dominates both metrics — S1=0.57 on peak T_in, S1=0.72 on % hours >26°C. The tunnel air temperature assumption is the single biggest uncertainty in the model.
- `D_CONC_EFF` is second on peak T_in (S1=0.28) but near-zero on comfort hours — thermal mass smooths peaks without changing how often the station runs hot.
- `AIRFLOW_OVERPRESSURE` is confirmed irrelevant on both metrics.

**Sobol B (water regime, 13 params, August week):**
- `T_CW_EXT_HYST` (chiller restart threshold) explains >90% of Q_water_cool variance — a pure regulation parameter choice.
- Water supply temperatures and glycol properties are noise.

---

## Status

**Week 3, Session 14.** Major model overhaul:
- Geometry: stair entry end wall excluded from A_soil. Concrete wall above PSD (UA_tun_wall = 192.5 W/K) added as separate tunnel-facing term.
- Physics: T_tun split into day (T_ext+10) / night (T_ext+5). Staircase passive infiltration (Q_stair = 1.25 m³/s) added to ODE.
- Setpoints: two-zone law replacing sliding law (heat 21°C / cool T_ext−6°C, dead band 15–20°C).
- Parameters: BASELINE_W 5000→500W (corridor only). T_COOL_FIXED replaced by T_COOL_DELTA=6°C.
- Sobol A/B/C rerun on updated model. All files cleaned — jargon moved to parameters.md, inline comments only in code.
- Next (new laptop): Sobol C rerun at N=1024, sweep on D_CONC_EFF × T_TUN_OFFSET_DAY, CO₂ emissions layer, psychrometric layer, Pareto front.