# Metro HVAC Optimizer

A physics-based digital twin of an underground metro platform in Paris, built to test and optimize HVAC regulation strategies on **public data only** — no proprietary files, no engineering-firm data.

It models one year of thermal and humidity dynamics on a real platform geometry, replicates a realistic AHU control strategy, ranks which control parameters actually matter (Sobol), finds the cost-vs-comfort sweet spot (Pareto), and serves the whole thing as an interactive web app.

**🔗 Live app:** https://metro-hvac--optimizer.streamlit.app/
**📦 Repo:** https://github.com/henrynasr/metro-hvac-optimizer

> The sections below explain the project end to end, physics included. If you just want the gist: it predicts a metro platform's energy, cost, CO₂ and comfort hour by hour, then tells you which knobs to turn. For the more curious, read on — every layer is broken down below, formulas and all.

---

## 1. The setting

Picture a side platform in a Paris underground metro stop. You walk down a staircase from street level onto a long, narrow concrete box. Trains arrive on one side behind platform screen doors (PSD).

That box is the zone we model. It measures **55 m long, 4 m wide and 4.2 m high** — a buried concrete corridor of roughly 924 m³. Each of its six faces behaves differently:

| Face | Boundary | Why |
|---|---|---|
| Back wall | **Soil** | Buried in stable ground |
| Ceiling | **Soil** | Buried |
| Floor | **Soil** | Buried |
| One end wall | **Soil** | Closed end |
| Tunnel side | **Tunnel air** | PSD glass facade + concrete wall above it |
| Stair side | **Open** | Staircase up to the street |

So **four faces touch soil**, **one faces the tunnel**, and **one is open to the outside** via the stairs. This split is the backbone of the whole thermal model — each boundary pushes or pulls heat at a different temperature.

---

## 2. How heat moves — the ODE

The platform air temperature `T_in` evolves over time following a **lumped-capacitance energy balance**: the whole zone is treated as one well-mixed thermal mass with capacitance `C`, and every heat flow in or out changes its temperature.

```
C · dT_in/dt = (UA_facade + UA_tun_wall + ρcp·V̇_inf)·(T_tun − T_in)   ← tunnel side
             + UA_soil·(T_soil − T_in)                                  ← soil sides
             + Q_int                                                    ← people + equipment
             + ρcp·Q̇_stair·(T_ext − T_in)                              ← fresh air via stairs
             + Q_curtain_zone                                           ← air curtain heat spill
             − Q_hvac                                                   ← what the AHU adds/removes
```

In plain words: the rate of temperature change equals the sum of all heat gains and losses, divided by the thermal mass. A bigger mass (`C`) means the temperature moves more slowly — the platform's concrete acts as a thermal flywheel, with a time constant of roughly **35 hours**.

Each term:

- **Tunnel side.** The facade (PSD glass) and the concrete wall above it both face the tunnel, not the outside. Tunnel air is hotter than outdoor air because of train braking, motor losses and passenger heat. We model it as `T_tun = T_ext + 10°C` during service, `+5°C` at night when trains stop. The `ρcp·V̇_inf` term is **door-opening infiltration**: every time a train stops and the PSD opens, a slug of tunnel air exchanges with the platform. Its size depends on train frequency (headway) and time of day.

- **Soil side.** The four buried faces exchange heat with ground at a near-constant `T_soil ≈ 15°C` — deep enough that seasonal swings are damped out. In winter the soil warms the platform; in summer it cools it.

- **Internal gains** `Q_int`. Sensible heat from people plus equipment (lighting, screens). More on the occupancy model in §5.

- **Staircase fresh air.** Outdoor air flows in passively through the open staircase. This is modulated (see §4) — fully open in mild weather, throttled by an air curtain when it's cold, shut behind a metal door at night.

- **Air curtain spill** `Q_curtain_zone`. When the curtain runs, part of its heat lands inside the platform (see §6).

- **`Q_hvac`.** The AHU's contribution — heating, cooling, or pure ventilation. Sign convention: positive `Q_hvac` = heat removed from the zone. This is where regulation lives (§3).

The ODE is integrated over a full year (8,784 hourly steps) with `scipy.integrate.solve_ivp`.

---

## 3. Regulation — the control strategy

A real AHU doesn't hold one fixed temperature. It follows an **outdoor-compensated sliding law**: the colder it is outside, the warmer the target inside, within limits. This avoids overheating in shoulder seasons and saves energy.

### Setpoint law (5 zones, by outdoor temperature)

```
T_ext ≤ 5°C        → heat to 18°C
5°C < T_ext < 15°C → heat to 18–20°C (linear ramp)
15°C ≤ T_ext ≤ 22°C → dead band — no heating, no cooling, pure ventilation
22°C < T_ext ≤ 32°C → cool to 26°C
T_ext > 32°C       → cool to (T_ext − 6°C)
```

The **dead band** is the key idea: between 15 and 22°C outside, the AHU spends nothing on heating or cooling — it just blows outdoor air through for air quality. At night (01h–05h) the law drops to anti-freeze only (5°C), since the station is closed.

### Airflow

Airflow scales with occupancy (more people → more fresh air for CO₂/hygiene), clamped between a minimum (overpressure + smoke-control floor) and a maximum. When heating or cooling load exceeds what the hygiene airflow can carry, the AHU **boosts** flow up to its max to meet the load.

### Water circuits + hysteresis

Heating and cooling are delivered by **water coils**, not direct electric:

- **Hot water** slides from 50°C supply (at −7°C outdoor) down to 40°C (at +15°C outdoor). The circuit **shuts off above 15°C** outdoor and **restarts at 13°C** — that gap is **hysteresis**, and it stops the circuit from rapidly switching on/off (short-cycling) when the outdoor temperature hovers near the threshold.
- **Cold water** is fixed: 8°C supply, 13°C return. It **shuts off below 26°C** outdoor and **restarts at 27°C** — same hysteresis logic.

---

## 4. Fresh air through the stairs

The open staircase is a free source of outdoor air — useful in summer, a liability in winter. It's modulated in three regimes:

| Condition | State | Effective flow |
|---|---|---|
| Mild weather (service hours) | Fully open | 100% |
| Cold weather (`T_ext` below threshold) | **Air curtain** active | ~35% |
| Night (01h–05h) | **Metal door** closed | ~8% (leakage) |

The air curtain is a jet of warm air across the stair opening that reduces cold-air infiltration without physically sealing the entrance.

---

## 5. Occupancy, internal gains and humidity sources

People are the dominant internal load, and they bring three things: **sensible heat** (warms the air), **latent heat** (moisture from breathing/sweat), and a need for **fresh air**.

Occupancy comes from **real RATP ridership data** (see §9), shaped into hourly profiles by day-type (working day / school-holiday weekday / weekend), anchored to a realistic peak headcount. Each person contributes a sensible heat load to `Q_int` and a moisture load to the humidity balance. Equipment adds a small fixed baseline.

---

## 6. The air curtain

A dedicated unit, modeled explicitly because it turned out to be a real chunk of the energy bill.

- It's a slot nozzle across the stair opening, blowing a warm jet downward.
- The jet has a defined cross-section, outlet velocity, and temperature rise above indoor air.
- **A fixed fraction of the curtain's heat stays inside the platform** (it spills into the zone — that's the `Q_curtain_zone` term in the ODE). The rest escapes up the stairwell and is lost.
- It runs on the same heat-pump hot-water branch as the main heating, so it shares the variable COP (§8).

---

## 7. Humidity and dehumidification

Temperature alone doesn't define comfort — humidity matters too. A separate **psychrometric layer** computes indoor relative humidity hour by hour.

The core quantity is the **humidity ratio** `W` (grams of water per kg of dry air). We track a moisture balance:

```
W_in = W_ext + (moisture from people) / (ρ · total fresh-air flow)
```

Indoor moisture is outdoor moisture plus what people add, diluted by all the fresh air coming in (stairs + infiltration + the outdoor fraction of AHU air). From `W_in` and `T_in` we back out relative humidity.

### Dehumidification by condensation

When cooling, air passes over a cold coil. If the **dew point of the mixed air exceeds the coil surface temperature**, water condenses on the coil — that's how the system dehumidifies. The condensed water carries away **latent heat**, which we add as an extra cooling load (and extra electricity).

### The winter problem

In winter, outdoor air is cold and dry. Once it's warmed to platform temperature, its relative humidity drops well below the comfort band — the air is simply too dry. **Condensation cooling can't fix dryness; only adding water can.** Without a dedicated humidifier (not installed in this model), winter dryness is a comfort floor we can't beat. This is stated honestly as a limitation, not engineered around.

---

## 8. Emissions, cost, and why the hour matters

Thermal loads become **electricity** through equipment efficiencies, and electricity becomes **CO₂** and **cost** — but neither conversion is constant through the day.

### Variable heating COP

Heating runs on an **air-source heat pump**. Its efficiency (COP) isn't fixed — it follows Carnot physics:

```
COP_heat = η_carnot · T_hw_K / (T_hw_K − T_ext_K)
```

The colder it is outside, the bigger the temperature lift the pump has to make, and the **worse** the COP. So winter heating is doubly expensive: more heat needed, delivered less efficiently. Cooling COP is treated as fixed (chilled water at a stable 8°C). Fan power follows the **affinity cube law** — `P_fan ∝ (flow)³` — so small airflow reductions yield large fan-energy savings.

### Time-varying CO₂

A kWh consumed at 3am isn't as dirty as one at 7pm. We use **RTE éCO2mix** data — the actual French grid carbon intensity (gCO₂/kWh), hour by hour, for 2024. Same energy, different CO₂ depending on when it's drawn.

### Time-varying price

Electricity is billed on a **horosaisonnier tariff** (a Tarif Jaune proxy): four price tiers split by season (winter/summer) and time of day (peak/off-peak). Off-peak winter nights are cheapest; peak winter days are most expensive.

---

## 9. Data sources (all public)

| Source | What | Use |
|---|---|---|
| **Open-Meteo** (ERA5 reanalysis) | 30 years of hourly Paris weather (temperature, humidity) | Outdoor boundary conditions — the CSV the model runs on |
| **RATP — Fréquentation du pôle La Défense** (IDFM open data) | Ridership counts by day-type and hour | Occupancy profiles — the Excel file driving the headcount |
| **RTE éCO2mix** | Hourly French grid CO₂ intensity, 2024 | Converting kWh → kgCO₂ at the right hour |

No proprietary, NDA, or engineering-firm data is used anywhere.

---

## 10. Comfort metrics

Comfort is scored on **service hours only** (the station is closed 01h–05h), in three tiers each for temperature and humidity:

| Tier | Temperature | Humidity |
|---|---|---|
| **Comfort** | 18–26°C | 40–60% RH |
| **Mild** | 14–18 or 26–28°C | 35–40 or 60–65% RH |
| **Discomfort** | <14 or >28°C | <35 or >65% RH |

**Combined** comfort = both temperature *and* humidity in their comfort band. **Combined discomfort** = at least one in its discomfort band. Everything in between is **mild**. The single headline number used for optimization is **% of service hours in combined discomfort** over the full year.

---

## 11. Sobol sensitivity analysis

With 15 controllable parameters, which ones actually drive the outcomes? A **Sobol global sensitivity analysis** answers this by variance decomposition: it attributes the variance in each output (cost, comfort, CO₂) to each input parameter.

- **S1 (first-order)** = how much that parameter matters on its own.
- **ST (total-order)** = its effect including interactions with other parameters.
- If S1 ≈ ST for everything and they sum to ~1, the system is **additive** (no interactions). A gap between S1 and ST signals interaction.

Saltelli sampling over the 15 parameters generates 16,384 full-year simulations. The ranked results — and how cost, comfort and CO₂ differ in what drives them — are **explorable interactively in the app** (Sobol page), so I won't freeze the numbers here.

---

## 12. Pareto front — optimizing cost vs comfort

You can't minimize cost and discomfort at the same time — cheaper usually means colder. The **Pareto front** is the set of configurations where you can't improve one without worsening the other.

Sweeping a grid of control levers (full-year sim for each), we extract the non-dominated configurations: the genuine trade-off frontier. Picking a point on it is a policy choice — how much comfort is worth how much money.

**Key insight: cost and CO₂ are tightly linked.** Both derive directly from energy consumption, so a configuration that's cheap is almost always low-carbon too — they rank nearly identically. The real tension is **energy vs comfort**, not cost vs CO₂. The Pareto explorer in the app lets you slide across configurations and see cost, CO₂ and the energy breakdown move together.

---

## 13. Pre-heating — a tested and rejected idea

Could we exploit the cheap night tariff by **banking heat in the concrete mass** overnight, to cut the expensive morning load? The thermal mass (τ ≈ 35h) is real, so it's worth checking. Two strategies were implemented and run for a full year:

1. **Setback override** — during winter nights (01h–05h), heat to 15°C instead of letting it drift to anti-freeze (5°C).
2. **Occupied HC boost** — during off-peak hours (22h–06h winter), overshoot the target to bank extra heat.

**Both lost on every metric** — more energy, more cost, no comfort gain. The reason: the zone naturally drifts to ~13–14°C overnight, and any target above that just spends energy that leaks out through staircase infiltration and ventilation before the morning peak arrives. The mass is a flywheel, but the losses are too large for it to work as a usable thermal battery. Chapter closed; the apparatus stays in the code as a documented dead end.

---

## 14. Fault detection (FDD) — catching a drifting sensor

Real sensors drift. If the outdoor temperature probe slowly biases high or low, every downstream decision degrades silently. This module detects that.

The idea: run the model on **clean** weather to get the "true" platform temperature (a stand-in for a trusted wall probe), then feed the controller a **corrupted** outdoor reading (here, +2°C bias starting July 1). Each hour, predict `T_in` from the suspect sensor and compare to truth. A **two-sided CUSUM** (cumulative sum) accumulates the residual and fires an alarm the instant it crosses a threshold.

The detector runs as a streaming loop — the exact same code would run unchanged on a live sensor feed. It catches the injected drift within a day, with no false alarms beforehand. It's **detection-only**: a single residual can't tell you *which* sensor drifted (that needs redundancy), and that's stated as future work.

---

## 15. The web app (Streamlit)

Everything above is served through an interactive **Streamlit** app, so it's not "clone the repo and run scripts" — it's a link you click.

**🔗 https://metro-hvac--optimizer.streamlit.app/**

Three pages:

1. **Live Simulation** — set the control parameters with sliders, hit run, get a full-year simulation (cost, energy, CO₂, discomfort + energy breakdown) in a few seconds. Real physics, computed on demand.
2. **Sobol Sensitivity** — the parameter rankings across cost, comfort and CO₂, switchable by metric.
3. **Pareto Explorer** — slide across precomputed configurations and watch the cost/comfort/CO₂ trade-off move against the optimal point.

---

## Project layout

| File | Role |
|---|---|
| `fetch_weather.py` | Pull 30y Paris weather from Open-Meteo |
| `constants.py` | Single source of truth — all parameters, units, Sobol ranges |
| `occupancy.py` | RATP profiles, day-type dispatch, infiltration rate |
| `regulation.py` | Setpoint law, staircase modulation, water circuits, the ODE |
| `simulation.py` | Full pipeline: ODE → HVAC → emissions → humidity → comfort |
| `emissions.py` | Electricity, CO₂, cost — variable COP, cube-law fans, hourly tariff |
| `humidity.py` | Psychrometric layer — RH, condensation, latent load |
| `sobol.py` | Sobol GSA — 15 params, 3 metrics |
| `pareto.py` | Pareto sweep — cost vs comfort |
| `compare.py` | Baseline vs optimized config, side by side |
| `preheat_comparison.py` | The rejected pre-heating experiment |
| `fdd.py` | Streaming CUSUM sensor-drift detector |
| `app.py` | The Streamlit web app |
| `utils.py` | Shared helpers |

### Setup

```bash
git clone https://github.com/henrynasr/metro-hvac-optimizer
cd metro-hvac-optimizer
python -m venv .venv
.venv\Scripts\activate          # Windows
pip install -r requirements.txt
python fetch_weather.py         # run once — pulls weather data
streamlit run app.py            # launch the app locally
```

---

## A note on assumptions

This is a **solo project on public data**. A production digital twin would be calibrated by engineering teams with measured values for every parameter. I didn't have that — so where a precise figure would have required on-site measurement or proprietary specs, I made a **documented, researched assumption** rather than fake precision.

The main ones: the **soil temperature** and the envelope **U-values** (how conductive the walls are), the **tunnel-side air temperature** offset (how much hotter than outdoor the tunnel runs), the **infiltration efficiency** per train stop, the **peak occupancy**, and the **fraction of air-curtain heat** retained inside the platform. Each is sourced from literature, standards, or comparable studies and flagged in `constants.py`. Getting them exact would need an instrumented station and a team — not worth it for a methodology demonstrator. The point here is the **method**, which transfers to any comparable platform once real values exist.