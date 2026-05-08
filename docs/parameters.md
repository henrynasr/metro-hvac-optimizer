# Energy Twin — Model Parameters

> Single reference document for all physical constants and modelling assumptions
> in the lumped-capacitance thermal model (one platform zone, one side, one AHU).
>
> **Scope:** 55 m × 4 m × 4.2 m platform zone. One of two symmetric zones.
> The two zones are not thermally coupled — model one, apply to both.
>
> **Conventions used throughout:**
> - **[SOBOL]** — parameter swept in Sobol sensitivity analysis. Range given.
> - **[ASSUMPTION]** — engineering estimate. No hard reference. Justified inline.
> - **[STUDY]** — rigorous computation requires a dedicated engineering study
>   (geotechnical, aeraulic, CFD, SES, structural thermal, etc.) that is
>   explicitly out of scope for a public-data model. Assumption taken in its
>   place and clearly flagged.
>
> *Session 11 — 2026-05-08. Public data only. No Fayat data.*

---

## 1. Geometry

| Parameter | Symbol | Value | Unit |
|---|---|---|---|
| Platform length | L | 55.0 | m |
| Platform width (one side) | W | 4.0 | m |
| Ceiling height | H | 4.2 | m |
| PSD facade height | H_f | 2.8 | m |
| Platform area | A_plat | 220.0 | m² |
| Platform air volume | V | 924.0 | m³ |
| Facade area (PSD wall) | A_f | 154.0 | m² |
| Buried envelope area (concrete/soil) | A_s | 704.6 | m² |

### Platform length — 55 m

GPE ligne 15 uses MP14 rolling stock in a 3-car formation.
Each car: ~18.5 m. Total train length: ~55.5 m → **55 m used**.

**Source:** SGP rolling stock technical note (MP14), referenced in the
Dossier d'Enquête Publique ligne 15 Sud (2015):
https://www.societedugrandparis.fr/gpe/ligne/ligne-15-sud

### Platform width — 4.0 m **[ASSUMPTION]**

Not independently measured. Derived from the occupancy target:
250 persons on 220 m² → 1.14 p/m² (Fruin LOS C).
220 m² / 55 m = **4.0 m**.
Width is a consequence of the occupancy assumption. No architectural drawing used.
**[SOBOL — range 3.5–5.0 m]**

### Ceiling height — 4.2 m **[ASSUMPTION]**

No GPE architectural drawing used. Typical for modern underground cut-and-cover
metro stations, consistent with published cross-sections of Paris M14 extensions
and London Jubilee Line deep-level stations.

A rigorous value requires the station architectural plans — a civil engineering
deliverable not available in public data for any GPE station.
**[SOBOL — range 3.8–5.0 m]**

### PSD facade height — 2.8 m **[ASSUMPTION]**

No GPE drawing used. Constrained by train door clearance (~2.1 m) plus structural
beam depth above. Standard for modern metro PSD assemblies.

Rigorous determination requires the PSD architectural and structural drawings —
a manufacturer and civil engineering deliverable not publicly available.
**[SOBOL — range 2.5–3.2 m]**

---

## 2. Thermal Envelope

The station is **fully buried**. Outdoor air temperature (T_ext) does NOT appear
directly as an envelope boundary condition. Two distinct thermal interfaces:

| Interface | Boundary temperature | Conductance |
|---|---|---|
| Facade (PSD glass wall) | T_tun (tunnel air) | UA_f = 1078 W/K (static) + dynamic infiltration |
| Structural box (concrete + soil) | T_soil = 15°C | UA_s = 289 W/K |

### 2a. Facade — U_f = 7.0 W/m²K

**Derivation (traceable chain):**

| Resistance component | Value (m²K/W) | Source |
|---|---|---|
| R_si — internal still air, platform side | 0.13 | ASHRAE Fundamentals 2017, Ch.26, Table 1 |
| R_glass — 12 mm tempered/laminated glass (λ = 1.0 W/mK) | 0.012 | EN ISO 10456:2007, Table 1 |
| R_se — tunnel side, moving air (piston effect) | 0.04 | ASHRAE Fundamentals 2017, Ch.26, Table 1 (forced convection) |
| **Total (glass only)** | **0.182** | → U_glass = 5.5 W/m²K |

Metal frame penalty (uninsulated aluminum/steel frames hold the glass panels):

- **ASHRAE 90.1-2019, Table A2.3:** "Uninsulated metal frame + single glazing"
  → assembly U = **6.0–8.0 W/m²K**
- **EN ISO 10077-1:2017** (thermal performance of windows, doors, shutters):
  equivalent range for uninsulated aluminum frames.

**U_f = 7.0 W/m²K** taken as mid-range, consistent with heavy-duty industrial
PSD specification.

A precise measured value requires a **certified laboratory thermal test** of the
actual PSD product under EN ISO 12567 (determination of thermal transmittance of
windows and doors). This is a manufacturer qualification test — not a public
document for any GPE PSD supplier.
**[SOBOL — range 5.5–8.5 W/m²K]**

### 2b. Structural walls — U_s = 0.41 W/m²K

**Derivation:**

| Resistance component | Value (m²K/W) | Source |
|---|---|---|
| R_si — internal still air | 0.13 | ASHRAE Fundamentals 2017, Ch.26 |
| R_concrete — 0.5 m wall, λ = 1.7 W/mK → R = 0.5/1.7 | 0.29 | EN ISO 10456:2007 (reinforced concrete) |
| R_soil — d_eff = 3.0 m, k_soil = 1.5 W/mK → R = 3.0/1.5 | 2.00 | ISO 13370:2017 §6; ASHRAE Fundamentals 2017, Ch.25 |
| **R_total** | **2.42** | → **U_s = 0.41 W/m²K** |

**Soil parameters:**
- k_soil = 1.5 W/mK: damp sandy clay, Paris basin geology, above water table.
  Consistent with values in ISO 13370 Annex A for temperate European urban soil.
- d_eff = 3.0 m: effective heat path through soil to the isothermal plane
  (undisturbed ground temperature). Value from ISO 13370 guidance for deep
  buried structures.

**R_soil sensitivity to soil moisture:**

| Condition | k_soil (W/mK) | R_soil (m²K/W) | U_s (W/m²K) |
|---|---|---|---|
| Dry soil | 1.0 | 3.0 | 0.35 |
| **Damp clay (used)** | **1.5** | **2.0** | **0.41** |
| Saturated / water table | 2.5 | 1.2 | 0.59 |

A rigorous U_s requires a **geotechnical site investigation**: soil boring logs,
in-situ thermal conductivity measurements (needle probe or TRT — Thermal Response
Test), and groundwater level monitoring over at least one annual cycle. This is a
dedicated geotechnical engineering study — not available from public data for any
specific GPE station site.
**[SOBOL — range 0.35–0.59 W/m²K]**

---

## 3. Boundary Temperatures

### T_soil = 15°C **[SOBOL — range 13–15°C]**

Undisturbed ground temperature at depth > 5 m, Île-de-France.
Measured range: **13–15°C**. Upper bound (15°C) used — conservative choice
that maximises the summer cooling load (worst-case design).

**Source:** BRGM (Bureau de Recherches Géologiques et Minières),
subsurface temperature atlas of metropolitan France:
https://www.brgm.fr/fr/enjeu/geothermie

### T_tun = T_ext + 5°C **[ASSUMPTION] [STUDY] [SOBOL — range 3–8°C, HIGH PRIORITY]**

Tunnel air is warmer than outdoor air due to heat dissipated by:
- Train braking (resistive energy, or regenerative losses when not recovered)
- Traction motor heat rejection
- Passenger body heat inside train cars

Range from the SES (Subway Environment Simulation) literature: **3–8°C**
above outdoor air, depending on headway, fleet thermal profile, and tunnel
ventilation effectiveness.

**No public measurement exists for GPE ligne 15 tunnels** (line not in commercial
operation at time of writing). The offset of 5°C is a mid-range engineering
estimate only.

A rigorous T_tun profile requires a **full SES aeraulic and thermodynamic
simulation** of the tunnel-station system: train speed profiles, motor heat
rejection curves (from traction system technical specifications), braking energy
balance, tunnel cross-section geometry, and tunnel ventilation capacity. This is
a dedicated engineering study typically conducted by the tunnel ventilation
specialist over 6–12 months during the detailed design phase. It is explicitly
out of scope for a public-data model.

This is the **highest-priority unknown** in the thermal boundary conditions.

---

## 4. Thermal Capacitance

| Component | Value | Unit |
|---|---|---|
| C_air (air mass only) | 1.115 × 10⁶ | J/K |
| C_concrete (effective 12 cm surface layer) | 1.718 × 10⁸ | J/K |
| **C_total** | **1.729 × 10⁸** | **J/K** |

Concrete dominates by **99%**. Air mass is negligible for a buried concrete box.

Static thermal time constant: τ = C_total / (UA_f + UA_s) = 1.729×10⁸ / 1367 ≈ **35 hours**.
The structural mass heats and cools very slowly. The AHU controls air temperature;
the concrete envelope drifts slowly around it.

**Air properties (standard):**
- ρ_air = 1.2 kg/m³, cp_air = 1005 J/(kg·K)
- Source: ASHRAE Fundamentals 2017, Chapter 1.

**Concrete properties:**
- ρ_conc = 2300 kg/m³ — EN 1991-1-1:2002
- cp_conc = 880 J/(kg·K) — EN ISO 10456:2007

### Effective concrete depth — d_eff = 0.12 m **[ASSUMPTION] [STUDY] [SOBOL — range 0.05–0.20 m]**

Only the surface layer of the concrete responds to short-term temperature changes
in the platform air. The thermal penetration depth is:

δ = √(α × t), where α = λ/(ρ×cp) = 1.7 / (2300 × 880) ≈ 8.4×10⁻⁷ m²/s

| Timescale | δ (m) |
|---|---|
| 1 h | 0.055 |
| 4 h | 0.110 |
| 24 h | 0.270 |

For diurnal station temperature cycles (1–4 h characteristic), **10–15 cm is
physically consistent**. Mid-range 12 cm used.

Reference convention: CIBSE Guide A (2015), Chapter 3 (thermal response and
dynamic simulation).

A rigorous effective depth requires either:
- **In-situ temperature sensors** at multiple depths in the structural concrete
  (a site instrumentation campaign), or
- A **2D/3D finite-element heat conduction model** of the station box (a
  dedicated structural thermal engineering study).

Both are out of scope here.

---

## 5. Occupancy

### Peak headcount — 250 persons/side

| | Value |
|---|---|
| Platform area | 220 m² |
| Peak density | 1.14 p/m² |
| LOS | C (busy, operationally acceptable) |

**Source:** Fruin, J.J. (1971). *Designing for Pedestrians: A Level-of-Service Concept*.
Port Authority of New York.

Referenced and tabulated in: **Metrolinx DS-12 Pedestrian Flow Modelling Design
Standard, v1.0, May 2024**, Section 5.2, Table 5-1.
URL: https://assets.metrolinx.com/image/upload/v1738255845/Documents/Engineering/DS-12_Pedestrian_Flow_Modelling_Design_Standard_v1.0.pdf

| LOS | m²/person | p/m² | Description |
|---|---|---|---|
| A | > 1.20 | < 0.83 | Free circulation |
| B | 0.90–1.20 | 0.83–1.10 | Restricted, comfortable |
| **C** | **0.65–0.90** | **1.10–1.54** | **Design target — busy, acceptable** |
| D | 0.28–0.65 | 1.54–3.57 | Crowded, short duration only |
| E | 0.19–0.28 | 3.57–5.26 | Near capacity |

**[SOBOL — range 150–350 persons (LOS B lower to LOS D lower)]**

### Sensible heat — 75 W/person

**Source:** ASHRAE Handbook of Fundamentals 2013, Chapter 18, Table 1.
*"Representative Rates at Which Heat and Moisture Are Given Off by Human Beings
in Different States of Activity"*

Activity: **"Standing, light work; walking"** at 24°C ambient → Sensible = 75 W, Total = 160 W.
Metro passengers (standing on platform, occasional slow walking) match this category.
**[SOBOL — range 70–90 W]**

### Latent heat — 60 W/person

Same source. Range 55–70 W at this activity level. Mid value used.
**Not in thermal ODE.** Reserved for the humidity submodel (future).

### Baseline equipment — 5000 W

LED platform lighting only. Escalators, gates, screens are concourse items — excluded.

**Derivation:**
- **Step 1 (T8 fluorescent reference):** SEAM4US project, Passeig de Gràcia Line 3
  station, Barcelona. Source: Ansuini et al. (2013), *Sustainable Energy Management
  for Underground Stations*, Portugal SB13, Table 1: 264 T8 lamps × 36 W = **9504 W**.
  URL: https://www.irbnet.de/daten/iconda/CIB_DC26414.pdf
- **Step 2 (LED correction):** Modern metro infrastructure uses LED throughout.
  LED delivers equivalent lux at 40–50% of T8 wattage. Factor: ×0.50.
- **Result:** 9504 × 0.50 = 4752 W → **5000 W**.

**[SOBOL — range 3000–8000 W]**

---

## 6. Infiltration

### Physical mechanism

Two components — both present even with full-height PSDs:

**A) Door-open bulk exchange (PSDs open, ~10 s per train stop)**
All 12 PSDs per side open simultaneously. The piston effect of the arriving
train has pressurised the tunnel ahead of it. On door opening, the transient
pressure differential drives a bi-directional air exchange — net infiltration
direction is tunnel→platform.

**B) Gap infiltration (PSDs closed)**
PSD assemblies have imperfect sealing: gaps at sliding door edges, top/bottom
frame joints, and service penetrations. Train piston pressure waves transiently
push tunnel air through these gaps. Platform overpressure (§7) keeps steady-state
flow direction outward, but cannot eliminate the transient spike at each arrival.

### Why the velocity method was rejected **[STUDY]**

The velocity method (V̇ = A_open × v_air × Δt) requires **v_air**: the air
velocity through the open door area during the piston pressure event. This cannot
be estimated without a **full aeraulic simulation of the tunnel-station interface**,
including train speed profile, tunnel cross-section geometry, PSD leakage
coefficient, and platform overpressure level. This is a dedicated SES
(Subway Environment Simulation) engineering study, conducted by the tunnel
ventilation specialist over several months during detailed design. It is
explicitly out of scope for a public-data model.

### Method used: Exchange Efficiency

V_cycle = η × V_platform
V̇_inf = V_cycle / T_headway

| Parameter | Symbol | Value | Basis |
|---|---|---|---|
| Exchange efficiency | η | 0.15 | Mid of 0.10–0.25 range |
| Volume per cycle | V_cycle | 138.6 m³ | η × 924 m³ |
| Peak V̇_inf | — | 1.155 m³/s | 138.6 / 120 s |
| Off-peak V̇_inf | — | 0.578 m³/s | 138.6 / 240 s |
| Night V̇_inf | — | 0.0 m³/s | No trains |

η = 0.15 **[ASSUMPTION] [STUDY]**

No GPE-specific or Paris Metro measurement found. Range 0.10–0.25 from PSD
subway measurement literature (search: *"Experimental study on air exchange
coefficient of subway stations with platform screen doors"*).

A rigorous η requires either:
- **Tracer gas measurements** during live train operations (SF₆ or CO₂ tracer,
  timed with train stops) — a dedicated measurement campaign requiring
  operational coordination with the metro operator; or
- A **full CFD/SES simulation** of the tunnel-station system with resolved
  piston flow and door leakage.

Both are dedicated engineering studies out of scope here.
**[SOBOL — range 0.10–0.25, HIGH PRIORITY — linear effect on Q_inf]**

### Train headways — GPE ligne 15 Ouest

**Source (peak):** Dossier d'Enquête d'Utilité Publique (DUP), ligne 15,
Pièce G — Caractéristiques de l'infrastructure:
https://www.societedugrandparis.fr/gpe/ligne/ligne-15-sud

| Regime | Headway | Trains/h | Basis |
|---|---|---|---|
| Peak (base model) | 2 min | 30 | DUP 2012 — committed value |
| Peak (Sobol upper) | 1.5 min | 40 | SGP 2022 ultimate capacity |
| Off-peak | 4 min | 15 | **[ASSUMPTION]** |
| Night | ∞ | 0 | No all-night service planned |

Off-peak and night headways: **[ASSUMPTION]** — the DUP gives peak headway only.
Derived from typical Île-de-France Mobilités metro scheduling practice.

### Night infiltration (residual gaps, no trains)

Treated as **zero** in this model. Real residual infiltration through PSD gaps
exists at all times but cannot be quantified without PSD leakage area data
(a manufacturer technical specification not publicly available).
Flagged for future refinement if gap infiltration is found to be significant
via Sobol analysis.

---

## 7. Overpressure — Minimum Platform Airflow

### Why overpressure is maintained even with full-height PSDs

Three independent physical requirements:

**1. Gap infiltration barrier**
PSD assemblies have imperfect sealing. Maintaining platform pressure above
tunnel pressure ensures net air flow direction is always platform→tunnel
through all gaps, preventing hot, humid, and particle-laden tunnel air
from entering through door edges, frame joints, and service penetrations.
If overpressure collapses, every gap becomes an infiltration point.

**2. Glass structural stress management**
A pressure differential across the PSD facade creates mechanical loading
on the glass panels and aluminum/steel frame connections. The overpressure
is sized to keep the differential within structural design limits — platform
slightly above tunnel, not the reverse. Reverse loading (tunnel > platform)
is the design failure mode (glass bending toward platform, door mechanism jamming).

**3. Smoke control (fire/emergency)**
In a fire event in the tunnel, combustion products must be directed away from
passengers. Platform overpressure is part of the smoke control strategy —
it ensures smoke stays tunnel-side and is extracted by dedicated tunnel fans,
not pushed onto the platform through PSD gaps.

### Minimum airflow — 2500 m³/h **[ASSUMPTION] [STUDY]**

The minimum airflow needed to sustain a target overpressure (typically 10–30 Pa
above tunnel) depends on:
- Total PSD leakage area (sum of all gap cross-sections across 12 door assemblies)
- Platform volume and geometry
- AHU pressure capability and duct layout

This requires a **full aeraulic pressure balance study** of the platform-PSD-tunnel
system, incorporating PSD leakage coefficients from the manufacturer and tunnel
pressure fluctuation data from the SES model. This is a dedicated ventilation
engineering deliverable — part of the SES study package and the fire safety
engineering report for the station. It is conducted by specialist ventilation
engineers and is not available from public data.

**2500 m³/h** is an engineering estimate consistent with typical practice for
PSD-equipped underground metro stations of this scale.
**[SOBOL — range 1500–4000 m³/h]**

---

## 8. Airflow Sizing

| Parameter | Value | Unit | Source |
|---|---|---|---|
| Regulatory minimum per person | 25 | m³/h/person | EN 16798-1:2019, Cat. II; ERP type GA |
| Overpressure minimum | 2500 | m³/h | Engineering estimate [ASSUMPTION] |
| Minimum airflow (zero occupancy) | 2500 | m³/h | Overpressure only |
| Maximum airflow (250 persons + margin) | 9625 | m³/h | (2500 + 250×25) × 1.10 |

**Additive logic (not max):**
The two requirements are physically independent — they cannot substitute for each
other. Overpressure is required at all times regardless of headcount; hygiene
demand scales with occupancy. Q_total = Q_overpressure + Q_occupancy.

---

## 9. Summary Table

| Parameter | Symbol | Value | Unit | Status |
|---|---|---|---|---|
| Platform length | L | 55.0 | m | Sourced (DUP/MP14) |
| Platform width | W | 4.0 | m | [ASSUMPTION] [SOBOL] |
| Ceiling height | H | 4.2 | m | [ASSUMPTION] [SOBOL] |
| Facade height | H_f | 2.8 | m | [ASSUMPTION] [SOBOL] |
| U_facade | U_f | 7.0 | W/m²K | Sourced (ASHRAE 90.1 / ISO 10077) [SOBOL] |
| U_soil | U_s | 0.41 | W/m²K | Derived (ISO 13370) [STUDY] [SOBOL] |
| T_soil | — | 15.0 | °C | Sourced (BRGM) [SOBOL] |
| T_tun offset | ΔT_tun | 5.0 | °C | [ASSUMPTION] [STUDY] [SOBOL HIGH] |
| C_total | C | 1.729×10⁸ | J/K | Derived (air + concrete) |
| Effective concrete depth | d_eff | 0.12 | m | [ASSUMPTION] [STUDY] [SOBOL] |
| Peak occupancy | N_peak | 250 | persons | Sourced (Fruin / Metrolinx DS-12) [SOBOL] |
| Sensible heat/person | q_s | 75 | W/person | Sourced (ASHRAE Fund. 2013 Ch.18) [SOBOL] |
| Baseline equipment | Q_base | 5000 | W | Sourced (SEAM4US + LED factor) [SOBOL] |
| Infiltration efficiency | η | 0.15 | — | [ASSUMPTION] [STUDY] [SOBOL HIGH] |
| Peak train headway | T_hw | 120 | s | Sourced (DUP 2012) |
| Overpressure airflow | Q_ovp | 2500 | m³/h | [ASSUMPTION] [STUDY] [SOBOL] |
| Regulatory air/person | q_vent | 25 | m³/h/p | Sourced (EN 16798-1) |

---

## 10. Parameters Not Yet Fixed

| Parameter | Status | Planned action |
|---|---|---|
| T_blow (AHU supply temperature) | Derived from setpoint + water regime | Output of dT_dt, not a fixed input |
| T_hot_water(T_ext) | Pending | Water regime law (S11) |
| T_cold_water(T_ext) | Pending | Water regime law (S11) |
| Q_hvac_max | Pending | Derived from water regime + airflow |
| PSD leakage area | Unknown (manufacturer spec) | Flagged — needed for rigorous Q_ovp |
| v_air during piston event | Unknown (SES study) | Flagged — needed for velocity method validation |

---

## 11. Sobol Sweep Priorities

Ordered by expected impact on model outputs (peak T_in, % hours > 26°C):

| Priority | Parameter | Range | Why |
|---|---|---|---|
| 1 | T_tun offset (ΔT_tun) | 3–8°C | Direct boundary condition on dominant facade term |
| 2 | η (infiltration efficiency) | 0.10–0.25 | Linear effect on Q_inf — large at peak headway |
| 3 | U_facade | 5.5–8.5 W/m²K | Multiplied by ΔT_tun — compounded uncertainty |
| 4 | U_soil / R_soil | 0.35–0.59 W/m²K | Lower priority: T_soil is stable |
| 5 | N_peak | 150–350 persons | Internal gains — already swept in S8 |
| 6 | d_eff (concrete depth) | 0.05–0.20 m | Affects C — expect low S1, moderate ST |
| 7 | T_soil | 13–15°C | Small range, quasi-steady — likely low sensitivity |

---

*Public data only. No Fayat data. No project-specific data.*
