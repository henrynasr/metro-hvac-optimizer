# Energy Twin — HVAC Sensitivity & Optimization Toolkit

Physics-based thermal model of a Paris underground metro stop — one side platform, one AHU.
Built on public data only. No proprietary files.

---

## What it does

Models the thermal dynamics of one buried platform zone using a lumped-capacitance ODE driven by real Paris weather (Open-Meteo ERA5) and real RATP occupancy data. The regulation layer replicates a real AHU control strategy: 5-zone setpoint law, airflow modulation with heating/cooling boost, water circuit scheduling, hysteresis logic, night setback, and staircase air curtain. A full-year humidity layer computes indoor RH, condensation risk, and latent cooling load.

A Sobol global sensitivity analysis (15 operator parameters, 16,384 runs) identifies the levers that matter. A Pareto front (6 levers, 1,296 configs) maps the cost-vs-comfort trade-off. A comparison script shows the before/after of the Pareto-optimal configuration.

**ODE:**

```
C · dT_in/dt = (UA_facade + UA_tun_wall + ρcp·V̇_inf)·(T_tun − T_in)
             +  UA_soil·(T_soil − T_in)
             +  Q_int
             +  ρcp·Q_stair(T_ext, hour)·(T_ext − T_in)
             +  Q_curtain_zone
             −  Q_hvac
```

- `T_tun = T_ext + 10°C` during service, `T_ext + 5°C` at night (01h–05h) — train braking, motor losses, passenger heat.
- `T_soil = 15°C` — stable ground temperature (BRGM). Applied to closed concrete side only.
- `UA_tun_wall` — concrete wall above PSD glass (1.4m × 55m), faces tunnel. U = 2.5 W/m²K.
- `V̇_inf` — dynamic PSD infiltration, f(hour, day-type), exchange efficiency method. Zero during night.
- `Q_stair(T_ext, hour)` — modulated staircase flow: full open (T_ext ≥ 7°C), curtain active (T_ext < 7°C, F=0.35), metal shutter night (F=0.08). Opening: 1.8m × 2.2m = 3.96 m².
- `Q_curtain_zone` — heat spill from air curtain into platform (40% of curtain coil power stays inside).
- `Q_int` — sensible heat from occupancy (75 W/person, ASHRAE) + equipment (500 W).
- `Q_hvac` — AHU power with airflow boost. T_blow: 15°C cooling, 30°C heating.

**Setpoint law (5-zone outdoor compensation):**

- T_ext ≤ 5°C → heat to 18°C (night: anti-freeze 5°C)
- 5 < T_ext < 15°C → heat to 18–20°C (linear ramp)
- 15 ≤ T_ext ≤ 22°C → dead band (pure ventilation)
- 22 < T_ext ≤ 32°C → cool to 26°C
- T_ext > 32°C → cool to T_ext − 6°C

**Variable COP (heating):** Carnot-based, air-source HP. `COP = η_carnot × T_hw_K / (T_hw_K − T_ext_K)`, η = 0.45, clamped [2.0, 8.0]. Same COP for air curtain (parallel branch on same HP). Cooling COP fixed (chilled water at 8°C).

---

## Data

- **Open-Meteo Historical Weather API** — ERA5 reanalysis, 30 years hourly Paris (1996–2025).
- **RATP — Fréquentation du pôle La Défense** (IDFM open data). Normalised hourly profiles by day-type (JOHV/JOVS/WKD), anchored to JOHV 18h = 250 persons.
- **RTE éCO2mix** — hourly grid CO₂ intensity (g/kWh), 2024 annual file.

---

## Scripts

| File | Role |
|---|---|
| `fetch_weather.py` | Pull 30y Paris weather from Open-Meteo |
| `constants.py` | Single source of truth — all parameters with units and Sobol ranges |
| `occupancy.py` | RATP profiles, day-type dispatch, `v_inf_m3s` |
| `regulation.py` | Setpoint law, staircase modulation, `dT_dt`, `build_Q_hvac_array`, water regime |
| `simulation.py` | Full pipeline: ODE → HVAC → emissions → humidity → latent correction → comfort classification → plots |
| `emissions.py` | Electricity, CO₂, cost — variable COP heating, cube law fans, RTE éCO2mix, horosaisonnier tariff |
| `humidity.py` | Psychrometric layer — RH_in, condensation risk, latent cooling load |
| `sobol.py` | Sobol GSA — 15 operator params, Saltelli N=512, 3 metrics (cost, comfort, CO₂) |
| `pareto.py` | Pareto sweep — 6 levers, 1,296 configs, cost vs combined discomfort |
| `preheat_compare.py` | Pre-heat strategy comparison — none / setback_override / hc, vs baseline |
| `compare.py` | Baseline vs Pareto-optimal — side-by-side comparison |
| `utils.py` | `load_data`, `style_axes` |
| `docs/parameters.md` | Full parameter table — values, sources, derivations |

---

## Setup

```bash
git clone https://github.com/henrynasr/energy-twin
cd energy-twin
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
python fetch_weather.py         # run once — pulls weather data
python simulation.py            # baseline annual run + plots (~2 min)
python compare.py               # baseline vs optimal (~4 min)
python pareto.py                # 1,296 configs (~43 min)
python sobol.py                 # 16,384 configs (~9 hours)
```

---

## Key findings

### Sobol GSA (15 operator parameters, N=512, 16,384 runs, full year 2024)

Three metrics: annual cost, combined discomfort (T+RH), annual CO₂.

**Cost and CO₂ are driven by the same physics** — nearly identical rankings, both additive (Σ S1 ≈ 0.96). `T_HEAT_LOW_C` alone explains 86% of cost variance.

**Combined discomfort is non-additive** (Σ S1 = 0.65). `T_HEAT_LOW_C` still dominates (S1=0.36, ST=0.60) but interacts with other parameters. Five other levers have ST > 0.05: `T_HW_EXT_HYST`, `AIRFLOW_PER_PERSON`, `T_DEAD_LOW`, `T_BLOW_HEAT`, `FRAC_RETURN_AIR`.

**Dead parameters:** `T_DEAD_HIGH_C`, `T_NIGHT_SETBACK_C`, `T_BLOW_COOL_C` — flat zero on all metrics.

### Pareto front (6 levers, 1,296 configs, 84 non-dominated)

Axes: annual cost (€) vs combined discomfort (% service hours with T or RH outside comfort).

**Three regimes, driven by `T_HEAT_LOW_C`:**
- H=14: 2,161–2,875€, 25–30% discomfort. Cheap but cold — 5–8% temperature discomfort.
- H=16: 2,986–3,318€, 17–18% discomfort. Temperature nearly solved (1.5%). Remaining discomfort is humidity.
- H=18: 3,922–4,373€, 12.6–14.3% discomfort. Expensive, diminishing returns.

**HW=40°C on every front point.** Lower supply = better COP = always optimal.

**T_HW_EXT_HYST has zero effect** — every front point appears in triplicate (11/13/15°C, identical results). Locked at baseline.

**Comfort floor is ~12% and is entirely humidity** — dry winter air (<35% RH) and humid summer air (>65% RH). No operator setpoint can fix this without a humidifier/dehumidifier.

### Baseline vs Pareto-optimal

Optimal config: H=16, A/pp=18, HW=40, Stair=7, Hyst=13, Blow=28.

| Metric | Baseline | Optimal | Delta |
|---|---|---|---|
| Annual cost | 4,520€ | 2,986€ | **-1,534€ (-34%)** |
| Annual energy | 26,590 kWh | 17,566 kWh | -9,024 kWh (-34%) |
| Annual CO₂ | 712 kg | 477 kg | -235 kg (-33%) |
| Heating energy | 15,010 kWh | 7,151 kWh | -7,858 kWh (-52%) |
| Combined comfort (both OK) | 53.1% | 27.2% | -25.9pp |
| Combined mild (tolerable) | 34.0% | 54.8% | +20.8pp |
| Combined discomfort (≥1 bad) | 12.9% | 18.0% | +5.1pp |

**Trade-off:** save 1,534€/year by letting the platform sit at 14–18°C (mild zone) instead of 18–26°C for more hours. Temperature discomfort barely changes (1.1% → 1.7%). The 5pp increase in combined discomfort is mostly humidity — cooler air holds less moisture in winter. For a metro platform where passengers spend under 5 minutes, the mild zone is acceptable.

### Comfort tiers (three-tier classification)

**Temperature:** comfort 18–26°C, mild 14–18°C or 26–28°C, discomfort <14°C or >28°C.
**Humidity:** comfort 40–60%, mild 35–40% or 60–65%, discomfort <35% or >65%.
**Combined:** comfort = both OK, discomfort = at least one bad, mild = everything else.

### Pre-heating under time-of-use pricing (no viable load-shift)

Tested whether banking heat in the concrete mass during cheap HC hours (winter night, 0.155 €/kWh) reduces the expensive morning HP load (0.215 €/kWh). Two strategies vs baseline (`preheat_compare.py`):

- **Night setback override** (01h–05h heated to 15 above anti-freeze): slightly worse on both energy and cost.
- **Occupied HC boost** (overshoot to 20°C during 22h–06h): **+1.1% energy, +0.6% cost** — strictly worse.

**Finding: no load-shifting opportunity exists on this geometry.** Any target above the zone's natural overnight drift (~13–14°C) costs more energy and more money — the heat is lost before the morning peak. The τ=35h thermal mass is real but the staircase infiltration and ventilation losses dominate it, so the mass is not a usable thermal battery. Pre-heating chapter closed; `PREHEAT_STRATEGY = "none"` retained as baseline. The strategy apparatus stays in `regulation.py` (`T_setpoint_preheat`) as the tested-and-rejected lever.

Note: switching from flat 0.17 €/kWh to horosaisonnier moved annual cost only +1.1% (4,520€ → 4,570€), confirming current consumption is not aligned with — nor exploitable by — the price signal.

---

## Known limitations

- Night window 23h–01h modeled as off-peak (4 min headway). Real last trains run ~15 min intervals.
- Air curtain 60% heat spill outside is lost, not recovered.
- Electricity tariff is horosaisonnier (Tarif Jaune proxy: HPH/HCH/HPE/HCE). HC window 22h–06h Mon–Sat assumed (Enedis sets it locally, not published generically).
- Coil contact factor = 1.0 (all air touches fins). Overestimates dehumidification ~10-20%.

---

## Status

**Week 6 — 2026-06-06.** Time-of-use pricing + pre-heating analysis.
- `emissions.py`: horosaisonnier tariff (Tarif Jaune proxy — HPH/HCH/HPE/HCE), `get_hourly_price()`, per-timestep cost. Flat→ToU moved annual cost only +1.1%.
- `regulation.py`: `T_setpoint_preheat()` wrapper, `PREHEAT_STRATEGY` lever (none / setback_override / hc). Pre-heat target during winter HC hours.
- `preheat_compare.py`: 3-strategy annual comparison via monkey-patch.
- **Finding: no viable pre-heating on this geometry.** Both strategies break even or lose (-0.1% to +1.1% energy). Infiltration-dominated zone — thermal mass not a usable battery. Chapter closed, `PREHEAT_STRATEGY = "none"` baseline.
- Next: TBD — options open (FDD stub, Monte Carlo on occupancy, SQL + Streamlit).

### Earlier
**Week 5 — 2026-05-30.** Sobol GSA + Pareto + comparison pipeline complete.
- `simulation.py` replaces `thermal_model.py`: full pipeline with humidity, latent correction, 3-tier comfort classification, 6 plot outputs.
- Sobol GSA: 15 operator params, N=512, 16,384 runs (8.6h). T_HEAT_LOW_C dominates cost (S1=0.86). Comfort is non-additive.
- Pareto: 6 levers, 1,296 configs, 84 non-dominated. Optimal: H=16, A=18, HW=40 → 2,986€/year.
- `compare.py`: baseline vs optimal side-by-side. -1,534€/year (-34%), +5.1pp discomfort.
- Deleted: `sweep.py` (superseded), `sobol_A.py`, `sobol_B.py` (replaced by unified `sobol.py`).