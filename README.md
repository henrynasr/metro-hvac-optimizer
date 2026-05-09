# Energy Twin — HVAC Sensitivity & Optimization Toolkit

Physics-based thermal model of an underground metro station platform (GPE Ligne 15 Ouest).
Built on public data only. No proprietary files.

---

## What it does

Models the thermal dynamics of one buried platform zone using a lumped-capacitance ODE driven by real Paris weather (Open-Meteo ERA5) and real RATP occupancy data. The regulation layer replicates a real AHU control strategy: setpoint law, airflow modulation, water circuit scheduling, and hysteresis logic. A Sobol global sensitivity analysis identifies the parameters that matter most.

**ODE:**

```
C · dT_in/dt = (UA_facade + ρcp·V̇_inf)·(T_tun − T_in)
             +  UA_soil·(T_soil − T_in)
             +  Q_int
             −  Q_hvac
```

- `T_tun = T_ext + 5°C` (tunnel air, piston effect + train heat)
- `T_soil = 15°C` (stable ground temperature, BRGM)
- `V̇_inf` = dynamic infiltration via PSD door openings, f(hour, day-type)
- `Q_int` = sensible heat from occupancy (75 W/person, ASHRAE) + equipment (5 kW LED)
- `Q_hvac` = AHU power. T_blow fixed: 15°C cooling, 30°C heating

**Water regime:**

Hot circuit: 50→35°C supply over T_ext −7→12°C. Off above 12°C (hysteresis: restarts at 10°C).
Cold circuit: 12→8°C supply over T_ext 26→31°C. Off below 26°C (hysteresis: restarts at 28°C).
Water flow derived from coil energy balance: `Q_water × ρ_glycol × Cp_glycol × ΔT_water = Q_air × ρcp_air × (T_blow − T_mix)`
where `T_mix = 0.7·T_in + 0.3·T_ext` (70% return air, industry practice).

---

## Data

- **Open-Meteo Historical Weather API** — ERA5 reanalysis, 30 years hourly Paris (1996–2025). Raw CSV gitignored.
- **RATP — Fréquentation du pôle La Défense** (IDFM open data). Normalised hourly profiles by day-type (JOHV/JOVS/WKD), anchored to JOHV 18h = 250 persons.

---

## Scripts

| File | Role |
|---|---|
| `fetch_weather.py` | Pull 30y Paris weather from Open-Meteo |
| `constants.py` | Single source of truth — all parameters, sourced or flagged |
| `occupancy.py` | RATP profiles, day-type dispatch, `v_inf_m3s` |
| `regulation.py` | `T_setpoint`, `airflow_total`, `dT_dt`, `build_Q_hvac_array`, water regime |
| `thermal_model.py` | Single-run simulation, 4×2 panel plot |
| `sweep.py` | 2D parameter sweep — pending update to new ODE |
| `sobol.py` | 6-parameter Sobol GSA via SALib — pending update to new ODE |
| `utils.py` | `load_data`, `style_axes` |
| `docs/parameters.md` | Full parameter table — sources, derivations, Sobol priority flags |

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

## Key findings so far

**July 2024:** T_in 21–25°C, T_ext 10–24°C. Dead-band ventilation (Q_vent 5–17 kW) dominates — tunnel air continuously pulls T_in down. Active cooling fires only on the hottest days. Water circuit correctly off (T_ext < 26°C throughout most of the week).

**Steady-state (no HVAC, mean July conditions):** T_eq ≈ 28.6°C — confirms active regulation is necessary in summer.

---

## Status

**Week 2, Session 12.** Water regime layer complete: `T_hot_water_supply(T_ext)` and `T_cold_water_supply(T_ext)` with 2°C hysteresis on both circuits. Q_water post-hoc output added to `build_Q_hvac_array`. T_mix (70/30 return/fresh air mix) integrated. Water circuit availability check added to both `dT_dt` and `build_Q_hvac_array` — HVAC power correctly zeroed when circuit is off. Next: update `sweep.py` and `sobol.py` to new split ODE, re-run Sobol on new parameter set (T_tun offset, η, U_facade, U_soil).
