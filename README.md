# Energy Twin — HVAC Sensitivity & Optimization Toolkit

Physics-based thermal model of a Paris underground metro stop — one side platform, one AHU.
Built on public data only. No proprietary files.

---

## What it does

Models the thermal dynamics of one buried platform zone using a lumped-capacitance ODE driven by real Paris weather (Open-Meteo ERA5) and real RATP occupancy data. The regulation layer replicates a real AHU control strategy: 5-zone setpoint law, airflow modulation with heating/cooling boost, water circuit scheduling, hysteresis logic, night setback, and staircase air curtain. A Sobol global sensitivity analysis identifies the parameters that matter most. A Pareto front over 5 regulation levers maps the energy-vs-comfort trade-off.

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
- `Q_hvac` — AHU power with airflow boost. T_blow fixed: 15°C cooling, 30°C heating. When water circuit is off but setpoint active, AHU blows at hygiene flow (pure ventilation, assigned to Q_vent).

**Setpoint law (5-zone):**

- T_ext ≤ 5°C → heat to 18°C (night 01h–05h: anti-freeze 5°C only)
- 5 < T_ext < 15°C → heat to 18–20°C (linear ramp)
- 15 ≤ T_ext ≤ 22°C → dead band (pure ventilation with outdoor air)
- 22 < T_ext ≤ 32°C → cool to 26°C
- T_ext > 32°C → cool to min(27, T_ext − 6)°C

**Water regime:**

Hot circuit: 50→40°C supply over T_ext −7→15°C. Off above 15°C (restarts at 13°C). ΔT = 5K constant.
Cold circuit: fixed 8°C supply, 13°C return, ΔT = 5K. Off below 26°C (restarts at 27°C).

**Variable COP (heating):**

Heat pump COP is Carnot-based, air-source: `COP = η_carnot × T_hw_K / (T_hw_K − T_ext_K)`, with `η_carnot = 0.45`. COP varies per timestep with hot water supply temperature and outdoor temperature. Clamped [2.0, 8.0]. Coldest hours (T_ext = -5°C, T_hw = 49°C): COP ≈ 2.9. Mild weather (T_ext = 12°C, T_hw = 41°C): COP ≈ 5.0. Same variable COP applies to air curtain (parallel branch on same heat pump). Cooling COP stays fixed (chilled water at constant 8°C).

**Air curtain:**

Dedicated hot-water coil above stair opening. Active during service when T_ext < 7°C. Nozzle: 0.08m × 2.5m = 0.2 m², jet at 7 m/s → 1.4 m³/s. Total electric ~10 kW (heater + fan). Fed by same heat pump as AHU heating circuit (COP = variable, see above). 40% of heat stays inside platform, 60% spills outside.

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
| `emissions.py` | Electricity, CO₂, cost — variable COP heating, cube law fans, RTE éCO2mix |
| `humidity.py` | Psychrometric layer — RH_in, condensation risk, latent cooling load |
| `pareto.py` | Multi-criteria Pareto sweep — 5 levers, 576 configs, front extraction + plot |
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
pip install psychrolib
python fetch_weather.py         # run once — pulls weather data
python thermal_model.py
python pareto.py                # ~23 min — 576 annual simulations
```

---

## Key findings

**Annual 2024 (one platform side, stair 1.8×2.2m, variable COP):**
T_in 13.7–29.1°C. E_total = 29,157 kWh, CO₂ = 779 kgCO₂, cost = 4,957€/year.
Heating = 17,616 kWh (60%), fans = 5,853 kWh (20%), curtain = 4,362 kWh (15%), cooling = 1,325 kWh (5%).
Comfort: 5.4% service hours above 26°C, 21.9% service hours below 18°C.

**Pareto front (576 configs, 5 levers, 37 non-dominated):**
- **HW=40°C on every front point.** Lower hot water supply = better COP = less electricity. 50°C water is never optimal.
- **Discomfort is almost entirely cold hours** — `%<18` ranges 2.7–51%, `%>26` stays 3–7%. The Pareto trade-off is about heating investment.
- **Sweet spot:** T_heat=21°C, T_cool=26°C, HW=40°C, overpressure=1500, airflow/person=25 → 31,838 kWh, 11.5% discomfort, 5,412€.
- **Comfort floor:** ~7.5% discomfort is the physical limit. Residual cold hours from staircase infiltration.
- **Steep improvement zone:** going from H=16 to H=20 buys 35pp of comfort for ~10,000 kWh. Beyond H=21, diminishing returns.

**Humidity (2024):**
RH comfort (40–60%): 64.0% of service hours. Too dry (<40%): 8.4% — winter, no fix without humidifier.
Too humid (>60%): 27.6% — dead band + heating hours, no active dehumidification.
Condensation risk: 135 hours. Latent cooling load: 3,014 kWh thermal.

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
- Electricity price fixed at 0.17 €/kWh. Real RATP contract is time-of-use.
- Coil contact factor = 1.0 (all air touches fins). Overestimates dehumidification ~10-20%.
- sobol_A.py and sobol_B.py not yet updated for S17 ODE signature changes.
- Pareto sweep uses monkey-patching (setattr on module namespaces) — functional but brittle. A parametric run_simulation() would be cleaner.

---

## Status

**Week 4, Session 21 — 2026-05-24.** Variable COP + Pareto front.
- Heating COP now Carnot-based, air-source. E_total +1.9% vs flat COP.
- Pareto front: 576 configs, 37 non-dominated. HW=40°C dominates. Sweet spot at 31,838 kWh / 11.5% discomfort.
- Next: full project review (physics + code), Sobol C rerun with humidity params, Pareto v2 (swap overpressure for T_STAIR_COLD_C), RH deep dive, cleanup.