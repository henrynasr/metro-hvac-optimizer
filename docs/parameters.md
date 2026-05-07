# Energy Twin — Model Parameters

> Parameter table for the lumped-capacitance thermal model of a single platform zone
> (one side, one AHU). All values validated or explicitly documented as engineering
> estimates during Session 10 (2026-05-08). Sobol sensitivity analysis will quantify
> the impact of uncertainty in each parameter — this table is the starting point, not
> the final word.
>
> Zone geometry: **4 m × 55 m × 7 m** (W × L × H), one platform side.
> Envelope area: 2×(55×7) + 2×(4×7) + 2×(4×55) = 770 + 56 + 440 = **1266 m²**
> Air volume: 4 × 55 × 7 = **1540 m³**

---

## 1. Occupancy — Peak headcount

| | Value |
|---|---|
| **Value used** | 250 persons/side |
| **Basis** | Fruin Level-of-Service framework |
| **Zone area** | 4 m × 55 m = 220 m² |
| **Resulting density** | 250 / 220 = **1.14 persons/m²** → LOS C (busy but acceptable) |

### Source

**Fruin, J.J. (1971). *Designing for Pedestrians: A Level-of-Service Concept*.**
Referenced and tabulated in: **Metrolinx DS-12 Pedestrian Flow Modelling Design Standard,
v1.0, May 2024**, Section 5.2 — Platforms & Queuing Areas.
URL: https://assets.metrolinx.com/image/upload/v1738255845/Documents/Engineering/DS-12_Pedestrian_Flow_Modelling_Design_Standard_v1.0.pdf

The Fruin LOS framework is the universal standard for metro platform density used
globally (TCRP Report 165, NFPA 130, Network Rail, London Underground). Density
values confirmed across multiple independent sources citing the same Fruin data.

### LOS reference table (waiting/queuing areas)

| LOS | m²/person | persons/m² | Description |
|---|---|---|---|
| A | > 1.2 | < 0.83 | Free circulation |
| B | 0.9–1.2 | 0.83–1.1 | Restricted but comfortable |
| **C** | **0.65–0.9** | **1.1–1.5** | **Design target — busy, acceptable** |
| D | 0.28–0.65 | 1.5–3.6 | Crowded, short durations only |
| E | 0.19–0.28 | 3.6–5.3 | Near capacity |

Our value of 1.14 p/m² sits at the LOS B/C boundary — a realistic peak design condition
for a mid-size GPE station, not an extreme scenario.

---

## 2. Envelope thermal conductance — UA

| | Value |
|---|---|
| **Value used** | UA = **600 W/K** |
| **U assumed** | 0.5 W/m²K |
| **A (envelope area)** | 1266 m² |
| **UA = U × A** | 0.5 × 1266 = 633 W/K → rounded to **600 W/K** |

### Source

**No published source gives a directly measured U-value for a buried metro station wall.**
An extensive literature search was conducted (Sadokierski & Thiffeault 2007 arXiv:0709.1748;
Shi et al. 2018 doi:10.1016/j.tust.2018.03.019; ASHRAE 90.1-2010 Table C6.10.1 via
EnergyPlus Engineering Reference https://bigladdersoftware.com/epx/docs/8-0/engineering-reference/page-026.html).
None applies directly to a fully buried station box.

**U = 0.5 W/m²K is an engineering estimate**, consistent with:
- Reinforced concrete wall (0.5 m, λ ≈ 1.7 W/mK) → R_concrete = 0.29 m²K/W
- Saturated Paris basin soil thermal resistance → R_soil significant but depth-dependent
  and not reliably sourced for this geometry
- Combined U = 1/(R_concrete + R_soil) → order of magnitude 0.3–0.8 W/m²K

**U = 0.5 W/m²K is taken as a central estimate pending calibration.**
Sobol analysis will sweep UA over ±50% to quantify sensitivity.

---

## 3. Thermal capacitance — C

| | Value |
|---|---|
| **Value used** | C = **1.86 × 10⁶ J/K** |
| **Represents** | Air mass in the platform zone only |
| **ρ_air** | 1.2 kg/m³ |
| **cp_air** | 1005 J/(kg·K) |
| **V_air** | 1540 m³ |
| **C = ρ × cp × V** | 1.2 × 1005 × 1540 = **1.86 × 10⁶ J/K** |

### Source

Standard dry air properties at ~20°C, atmospheric pressure. These are fundamental
thermodynamic constants, not model assumptions:
- ρ_air = 1.2 kg/m³: consistent with ideal gas law at 20°C, 1 atm
- cp_air = 1005 J/(kg·K): standard value, tabulated in any engineering thermodynamics
  reference (e.g. ASHRAE Fundamentals 2017, Chapter 1)

### Modelling note

T_in in the ODE represents **air temperature only**. Concrete thermal mass is not
included in C — the station structure is treated as a boundary condition (via UA),
not a lumped capacitance. This simplification means the model has low thermal
inertia (τ = C/UA = 1.86×10⁶/600 ≈ **3100 s ≈ 52 min**). The air heats and cools
quickly; the regulation layer provides the smoothing. This choice will be tested
via Sobol — a 2-node RC model (air + structure) is a natural next iteration if
sensitivity to C is found to be high.

---

## 4. Sensible heat per occupant

| | Value |
|---|---|
| **Value used** | **75 W/person (sensible only)** |
| **Activity** | Standing, light work / walking slowly |
| **Room temperature** | 24°C (ASHRAE reference condition) |
| **Latent (not in ODE)** | 55–70 W/person → goes into humidity model |

### Source

**ASHRAE Handbook of Fundamentals 2013, Chapter 18, Table 1**
*"Representative Rates at Which Heat and Moisture Are Given Off by Human Beings
in Different States of Activity"*

Confirmed from Scribd scan of the actual table:
https://www.scribd.com/document/686101040/05-ATTACHMENT-8-3-ASHRAE-Fundamentals-2013-Heat-Gain-from-Occupants

Relevant table entries (at 24°C room temperature, adult male):

| Activity | Location | Total W | Sensible W | Latent W |
|---|---|---|---|---|
| Standing, light work; walking | Department store | 160 | 75 | 55 (adj.) |
| Walking, standing | Drug store, bank | 160 | 75 | 70 |

Metro passengers (standing on platform, occasional slow walking) match the
"standing, light work; walking" category. At 26–27°C (warmer station), ASHRAE
notes sensible drops ~20% — value of 75 W is therefore slightly conservative,
which is appropriate for a design model.

### Modelling note

Only sensible heat (75 W/person) enters the thermal ODE as Q_internal.
Latent heat raises air humidity but not air temperature — it will feed the
humidity submodel when built.

---

## 5. Baseline equipment load

| | Value |
|---|---|
| **Value used** | **5 kW** (platform zone only) |
| **Covers** | LED platform lighting |
| **Excludes** | Escalators, ticket gates, screens (concourse zone, not platform) |

### Source and derivation

**Step 1 — T8 fluorescent baseline:**
SEAM4US project, Passeig de Gràcia Line 3 station, Barcelona (2013).
Source: Ansuini et al., *Sustainable Energy Management for Underground Stations:
Lighting Upgrade*, Portugal SB13 Conference Proceedings, pp. 347–354.
URL: https://www.irbnet.de/daten/iconda/CIB_DC26414.pdf
Fetched and verified. Table 1 of the paper (p.349 of proceedings):
**Platform: 264 T8 lamps × 36 W = 9,504 W ≈ 9.5 kW**

**Step 2 — LED correction:**
GPE stations use LED lighting throughout (standard for all new European metro
infrastructure post-2015). LED produces equivalent lux output at 40–50% of T8
wattage. This is confirmed directionally by the WMATA LED upgrade program
(48 stations, 17,000,000 kWh/year saved):
https://www.mcdean.com/making-dcs-metro-brighter-safer-and-more-sustainable/

**Step 3 — Result:**
9.5 kW × 0.50 = **~4.8 kW → rounded to 5 kW**

This is a platform-zone estimate only. A full station model would add:
escalators (~8–15 kW each, concourse), ticket gates (~150 W each), info screens
(~75 W each), and service room loads.

---

## 6. Summary table

| Parameter | Symbol | Value used | Unit | Basis |
|---|---|---|---|---|
| Peak occupancy | N_peak | 250 | persons/side | Fruin LOS C, 220 m² platform zone |
| Envelope conductance | UA | 600 | W/K | U=0.5 W/m²K (estimate) × A=1266 m² |
| Thermal capacitance | C | 1.86 × 10⁶ | J/K | Air mass only: ρ·cp·V |
| Sensible heat/person | q_sens | 75 | W/person | ASHRAE Fundamentals 2013, Ch.18, Table 1 |
| Latent heat/person | q_lat | 60 | W/person | ASHRAE Fundamentals 2013, Ch.18, Table 1 (humidity model only) |
| Baseline equipment | Q_base | 5,000 | W | SEAM4US Barcelona × LED correction factor |

---

## 7. Parameters NOT fixed here

| Parameter | Status | Next step |
|---|---|---|
| T_blow (supply air temp) | Derived, not fixed | Output of controller: T_blow = f(T_set, Q_air, T_in) |
| T_hot_water(T_ext) | Pending | S10/S11: water regime law with hysteresis |
| T_cold_water(T_ext) | Pending | S10/S11: water regime law with hysteresis |
| Q_hvac_max | Pending | Derived from water regime + airflow |

---

## 8. What comes next — Sobol on these parameters

The values above are starting points, not ground truth. The planned Sobol sensitivity
analysis (briefing question 1) will sweep each parameter over a defensible range and
quantify its contribution to variance in:
- Peak T_in
- % hours T_in > 26°C
- Annual HVAC energy

Parameters with low Sobol indices can be fixed confidently. Parameters with high
indices are the ones worth calibrating carefully against real data.

---

*Session 10 — 2026-05-08. Parameters validated one by one with source verification.
Public data only. No Fayat data.*
