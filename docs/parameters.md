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
> *Session 12 — 2026-05-09. Public data only. No Fayat data.*

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
**[SOBOL — range 3.8–5.0 m]**

### PSD facade height — 2.8 m **[ASSUMPTION]**

No GPE drawing used. Constrained by train door clearance (~2.1 m) plus structural
beam depth above. Standard for modern metro PSD assemblies.
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
| R_se — tunnel side, moving air (piston effect) | 0.04 | ASHRAE Fundamentals 2017, Ch.26, Table 1 |
| **Total (glass only)** | **0.182** | → U_glass = 5.5 W/m²K |

Metal frame penalty: ASHRAE 90.1-2019 Table A2.3 + EN ISO 10077-1:2017 →
"Uninsulated metal frame + single glazing" → assembly U = **6.0–8.0 W/m²K**.
**U_f = 7.0 W/m²K** taken as mid-range.
**[SOBOL — range 5.5–8.5 W/m²K]**

### 2b. Structural walls — U_s = 0.41 W/m²K

**Derivation:**

| Resistance component | Value (m²K/W) | Source |
|---|---|---|
| R_si — internal still air | 0.13 | ASHRAE Fundamentals 2017, Ch.26 |
| R_concrete — 0.5 m wall, λ = 1.7 W/mK | 0.29 | EN ISO 10456:2007 |
| R_soil — d_eff = 3.0 m, k_soil = 1.5 W/mK | 2.00 | ISO 13370:2017 §6 |
| **R_total** | **2.42** | → **U_s = 0.41 W/m²K** |

**R_soil sensitivity:**

| Condition | k_soil (W/mK) | U_s (W/m²K) |
|---|---|---|
| Dry soil | 1.0 | 0.35 |
| **Damp clay (used)** | **1.5** | **0.41** |
| Saturated | 2.5 | 0.59 |

**[SOBOL — range 0.35–0.59 W/m²K]**

---

## 3. Boundary Temperatures

### T_soil = 15°C **[SOBOL — range 13–15°C]**

Undisturbed ground temperature at depth > 5 m, Île-de-France. Upper bound used — conservative (maximises summer cooling load).

**Source:** BRGM subsurface temperature atlas: https://www.brgm.fr/fr/enjeu/geothermie

### T_tun = T_ext + 5°C **[ASSUMPTION] [STUDY] [SOBOL — range 3–8°C, HIGH PRIORITY]**

Tunnel air warmer than outdoor due to train braking losses, traction motor heat, passenger body heat.
No public measurement for GPE ligne 15 tunnels. Rigorous value requires full SES aeraulic and thermodynamic simulation — out of scope.

---

## 4. Thermal Capacitance

| Component | Value | Unit |
|---|---|---|
| C_air | 1.115 × 10⁶ | J/K |
| C_concrete (12 cm effective surface layer) | 1.718 × 10⁸ | J/K |
| **C_total** | **1.729 × 10⁸** | **J/K** |

Concrete dominates 99%. τ_static = C_total / (UA_f + UA_s) ≈ **35 hours**.

### Effective concrete depth — 0.12 m **[ASSUMPTION] [STUDY] [SOBOL — range 0.05–0.20 m]**

Thermal penetration depth δ = √(α·t), α = 8.4×10⁻⁷ m²/s:

| Timescale | δ (m) |
|---|---|
| 1 h | 0.055 |
| 4 h | 0.110 |
| 24 h | 0.270 |

10–15 cm consistent with diurnal cycles. 12 cm used (mid-range).
Reference: CIBSE Guide A (2015), Chapter 3.

---

## 5. Occupancy

### Peak headcount — 250 persons/side

Fruin LOS C at 1.14 p/m² on 220 m².

**Source:** Fruin (1971) via Metrolinx DS-12 v1.0, May 2024, §5.2 Table 5-1.
URL: https://assets.metrolinx.com/image/upload/v1738255845/Documents/Engineering/DS-12_Pedestrian_Flow_Modelling_Design_Standard_v1.0.pdf

**[SOBOL — range 150–350 persons]**

### Sensible heat — 75 W/person

**Source:** ASHRAE Fundamentals 2013, Ch.18, Table 1. "Standing, light work; walking" at 24°C.
**[SOBOL — range 70–90 W]**

### Latent heat — 60 W/person

Same source. Not in thermal ODE — reserved for psychrometric model (future).

### Baseline equipment — 5000 W

LED platform lighting only. Escalators and gates are concourse items — excluded.

Derivation: SEAM4US Barcelona (Ansuini et al. 2013), Table 1: 264 T8 × 36 W = 9504 W → ×0.50 LED factor = **5000 W**.
URL: https://www.irbnet.de/daten/iconda/CIB_DC26414.pdf
**[SOBOL — range 3000–8000 W]**

---

## 6. Infiltration

Exchange efficiency method: V̇_inf = η × V_platform / T_headway.

η = 0.15 **[ASSUMPTION] [STUDY] [SOBOL — range 0.10–0.25, HIGH PRIORITY]**

No GPE or Paris Metro measurement found. Range from PSD subway literature.
Rigorous η requires tracer gas measurement campaign or full CFD/SES simulation — out of scope.

| Regime | Headway | V̇_inf (m³/s) |
|---|---|---|
| Peak | 120 s (DUP 2012) | 1.155 |
| Off-peak | 240 s [ASSUMPTION] | 0.578 |
| Night | ∞ | 0.0 |

---

## 7. Overpressure — Minimum Platform Airflow

Maintained for: (1) gap infiltration barrier, (2) glass structural stress management, (3) smoke control.

**2500 m³/h per zone [ASSUMPTION] [STUDY]** — engineering estimate. Rigorous value requires aeraulic pressure balance study with PSD leakage coefficients — out of scope.
**[SOBOL — range 1500–4000 m³/h]**

---

## 8. Airflow Sizing

| Parameter | Value | Unit | Source |
|---|---|---|---|
| Regulatory minimum per person | 25 | m³/h/person | EN 16798-1:2019, Cat. II |
| Overpressure minimum | 2500 | m³/h | [ASSUMPTION] |
| AIRFLOW_MIN | 2500 | m³/h | Overpressure only (zero occupancy) |
| AIRFLOW_MAX | 9625 | m³/h | (2500 + 250×25) × 1.10 |

Additive logic: Q_total = Q_overpressure + Q_occupancy. Not max() — the two requirements are independent.

---

## 9. Water Regime

### Hot water circuit

Supply temperature: 50°C at T_ext = −7°C → 35°C at T_ext = 12°C (linear, `np.interp`).
Shut off above 12°C. Restarts when T_ext drops back to 10°C (2°C hysteresis band).

**Source (setpoints):** Fayat Energies Services internship documentation (Henry Nasr, 2026).
**Source (hysteresis):** ASHRAE Handbook — HVAC Systems and Equipment (2020), Ch.42 — deadband control in hydronic heating/cooling plant.

### Cold water circuit

Supply temperature: 12°C at T_ext = 26°C → 8°C at T_ext = 31°C (linear, `np.interp`).
Shut off below 26°C. Restarts when T_ext rises back to 28°C (2°C hysteresis band).
Same source for hysteresis logic.

### AHU coil energy balance

```
Q_water [m³/s] = Q_air × ρcp_air × dT_air / (ρ_glycol × Cp_glycol × ΔT_water)
```

Where:
- `dT_air = T_blow − T_mix` (heating) or `T_mix − T_blow` (cooling)
- `T_mix = 0.7 × T_in + 0.3 × T_ext` (70% return air / 30% fresh air)
- `ΔT_water_heat = 5 K` (50/45°C or 35/30°C — constant across operating range)
- `ΔT_water_cool = 6 K` (8/14°C or 12/18°C — constant across operating range)

Q_water is a **post-hoc secondary output** — does not modify Q_hvac or dT_dt.

### Glycol properties (30% propylene glycol mix)

| Property | Value | Unit | Status |
|---|---|---|---|
| Cp_glycol | 3800 | J/(kg·K) | [ASSUMPTION — standard value] |
| ρ_glycol | 1045 | kg/m³ | [ASSUMPTION — standard value] |

---

## 10. AHU Air Mix

**FRAC_RETURN_AIR = 0.70** (70% return air / 30% outdoor air)

T_mix = 0.70 × T_in + 0.30 × T_ext

No international standard (ASHRAE 62.1, EN 16798-3, NFPA 130) prescribes a numeric cap on RA fraction; OA volume is set by occupancy and IAQ requirements. Platform air classified as ASHRAE 62.1 Class 2/3 (elevated PM, metallic aerosols) → recirculation within zone acceptable.

**Sources for 70% RA:** Seoul SCAP system (PMC 2022, doi:10.3390/ijerph192013302); Delhi Metro coach design (~73% RA, Indian Express 2020).

---

## 11. Humidity Targets

| Threshold | Value | Source |
|---|---|---|
| Comfort lower bound | 40% RH | ASHRAE 55; EN ISO 7730 |
| Comfort upper bound | 60% RH | ASHRAE 55; Delhi Metro ECS (ICTRAM 2018) |
| Low alert (dryness) | 30% RH | ASHRAE 55 lower bound |
| High alert (condensation risk) | 70% RH | UCL underground moisture study, Wei et al., BSE 2021 |
| Structural limit (mold) | 80% surface RH | ASHRAE 160, 30-day running average |

Dehumidification load not yet modelled — psychrometric layer deferred.
When outdoor air exceeds 60% RH, the AHU cooling coil must extract latent heat in addition to sensible. This will be addressed in the psychrometric extension.

---

## 12. Sobol Sweep Priorities

| Priority | Parameter | Range | Why |
|---|---|---|---|
| 1 | T_tun offset | 3–8°C | Direct boundary condition on dominant facade term |
| 2 | η (infiltration) | 0.10–0.25 | Linear effect on Q_inf at peak headway |
| 3 | U_facade | 5.5–8.5 W/m²K | Multiplied by ΔT_tun — compounded uncertainty |
| 4 | U_soil | 0.35–0.59 W/m²K | Lower priority: T_soil is stable |
| 5 | N_peak | 150–350 persons | Internal gains |
| 6 | d_eff | 0.05–0.20 m | Affects C — expect low S1 |
| 7 | T_soil | 13–15°C | Small range — likely low sensitivity |

---

*Public data only.*
