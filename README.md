# Energy Twin — HVAC Sensitivity & Optimization Toolkit

Physics-based thermal model of a Paris underground metro stop — one side platform, one AHU.
Built on public data only. No proprietary files.

---

## What it does

Models the thermal dynamics of one buried platform zone using a lumped-capacitance ODE driven by real Paris weather (Open-Meteo ERA5) and real RATP occupancy data. The regulation layer replicates a real AHU control strategy: 5-zone setpoint law, airflow modulation with heating/cooling boost, water circuit scheduling, hysteresis logic, night setback, and staircase air curtain. A Sobol global sensitivity analysis identifies the parameters that matter most.

**ODE:**

```
C · dT_in/dt = (UA_facade + UA_tun_wall + ρcp·V̇_inf)·(T_tun − T_in)
             +  UA_soil·(T_soil − T_in)
             +  Q_int
             +  ρcp·Q_stair(T_ext, hour)·(T_ext − T_in)
             +  Q_curtain_zone
             −  Q_hvac
```

- `T_tun = T_ext + 10°C` during service, `T_ext + 5°C` at night (01h–05h) — train braking, motor losses, passenger heat. [SOBOL 5–15 / 3–8]
- `T_soil = 15°C` — stable ground temperature (BRGM). Applied to closed concrete side only.
- `UA_tun_wall` — concrete wall above PSD glass (1.4m × 55m), faces tunnel. U = 2.5 W/m²K.
- `V̇_inf` — dynamic PSD infiltration, f(hour, day-type), exchange efficiency method. Zero during night (01h–05h, no trains).
- `Q_stair(T_ext, hour)` — modulated staircase flow: full open (T_ext ≥ 7°C), curtain active (T_ext < 7°C, F=0.35), metal shutter night (F=0.08). Opening: 1.8m × 2.2m = 3.96 m², velocity 0.5 m/s → Q_base = 1.98 m³/s.
- `Q_curtain_zone` — heat spill from air curtain into platform (40% of curtain coil power stays inside).
- `Q_int` — sensible heat from occupancy (75 W/person, ASHRAE) + equipment (500 W LED/screens, corridor only).
- `Q_hvac` — AHU power with airflow boost. T_blow fixed: 15°C cooling, 30°C heating.

**Setpoint law (5-zone):**

- T_ext ≤ 5°C → heat to 18°C (night 01h–05h: anti-freeze 5°C only)
- 5 < T_ext < 15°C → heat to 18–20°C (linear ramp)
- 15 ≤ T_ext ≤ 22°C → dead band (pure ventilation with outdoor air)
- 22 < T_ext ≤ 32°C → cool to 26°C
- T_ext > 32°C → cool to min(27, T_ext − 6)°C

**Water regime:**

Hot circuit: 50→40°C supply over T_ext −7→15°C. Off above 15°C (restarts at 13°C). ΔT = 5K constant.
Cold circuit: fixed 8°C supply, 13°C return, ΔT = 5K. Off below 26°C (restarts at 27°C).

**Air curtain:**

Dedicated electric unit above stair opening. Active during service when T_ext < 7°C. Nozzle: 0.08m × 2.5m = 0.2 m², jet at 7 m/s → 1.4 m³/s. Total electric ~10 kW (heater + fan, COP = 1.0). 40% of heat stays inside platform, 60% spills outside.

---

## Data

- **Open-Meteo Historical Weather API** — ERA5 reanalysis, 30 years hourly Paris (1996–2025). Raw CSV gitignored.
- **RATP — Fréquentation du pôle La Défense** (IDFM open data). Normalised hourly profiles by day-type (JOHV/JOVS/WKD), anchored to JOHV 18h = 250 persons.
- **RTE éCO2mix** — real-time grid CO₂ intensity (g/kWh), 2024 annual file.

---

## Scripts

| File | Role |
|---|---|
| `fetch_weather.py` | Pull 30y Paris weather from Open-Meteo |
| `constants.py` | Single source of truth — all parameters with units and Sobol ranges |
| `occupancy.py` | RATP profiles, day-type dispatch, `v_inf_m3s` |
| `regulation.py` | `T_setpoint`, `q_stair_m3s`, `airflow_total`, `dT_dt`, `build_Q_hvac_array`, water regime |
| `thermal_model.py` | Single-run simulation, 4×2 panel plot |
| `emissions.py` | Electricity, CO₂, cost — post-hoc from Q arrays + RTE éCO2mix + curtain energy |
| `sobol_A.py` | Sobol A (screen) + C (5 survivors, N=512) via SALib. Self-contained parametric ODE. |
| `sobol_B.py` | Sobol B — water regime sensitivity, post-hoc. August week. |
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

**Annual 2024 (one platform side, stair 1.8×2.2m):**
T_in 15.4–29.2°C. E_total = 27,639 kWh, CO₂ = 780 kgCO₂, cost = 4,699€/year.
Heating = 13,138 kWh (48%), curtain = 10,087 kWh (36%), fans = 3,649 kWh (13%), cooling = 766 kWh (3%).
Comfort: 5.6% hours above 26°C, 18.1% hours below 18°C (service + night combined).

**Key finding — air curtain cost:** curtain runs 1,331 h/year and accounts for 36% of total electricity. At COP=1 (electric resistance), it is the single most expensive component after heating. Worth revisiting: dedicated heat pump unit or hot water feed from district heating would cut this significantly.

**Sobol C (5 params, N=512, July week):**
- `T_TUN_OFFSET_DAY` dominates both metrics — S1=0.57 on peak T_in, S1=0.72 on % hours >26°C.
- `D_CONC_EFF` is second on peak T_in (S1=0.28) but near-zero on comfort hours.
- `AIRFLOW_OVERPRESSURE` confirmed irrelevant.

**Sobol B (water regime, 13 params, August week):**
- `T_CW_EXT_HYST` (chiller restart threshold) explains >90% of Q_water_cool variance.
- Water supply temperatures and glycol properties are noise.

---

## Known limitations

- Night window 23h–01h modeled as off-peak (4 min headway). Real last trains run ~15 min intervals — V_inf overestimated in that window.
- Air curtain energy does not account for the 60% heat spill outside (lost, not recovered).
- COP fixed year-round. Real COP varies with supply temperature and part-load.
- Electricity price fixed at 0.17 €/kWh. Real RATP contract is time-of-use.
- sobol_A.py and sobol_B.py not yet updated for S17 ODE signature changes.

---

## Status

**Week 3, Session 17 — 2026-05-18.** Regulation overhaul + staircase modulation + air curtain.
- Dead band: was Q_hvac=0, now pure ventilation with outdoor air (physically correct).
- Heating + cooling airflow boost: AHU ramps up to AIRFLOW_MAX when hygiene flow can't meet load.
- Night setback: 01h–05h → anti-freeze 5°C only. Trains stop at 01h.
- Staircase: A_stair corrected to 1.8×2.2m = 3.96 m², Q_base = 1.98 m³/s. Three regimes: full open / curtain active (F=0.35) / metal shutter night (F=0.08).
- Air curtain: dedicated electric unit, ~10 kW, active when T_ext < 7°C during service. Heat spill (40%) added to ODE.
- Water: cold fixed 8/13°C (ΔT=5K). Hot min supply corrected to 40°C (return 35°C, ΔT=5K constant).
- emissions.py: curtain electric power (heater + fan, COP=1) added to all energy/CO₂/cost totals.
- Next: new laptop env setup, sobol_A/B update for new ODE signature, psychrometric layer, Pareto front. Revisit curtain energy source (COP=1 is expensive).