# =============================================================================
# constants.py — Energy Twin model parameters
# Single source of truth. Import everywhere with: from constants import *
#
# Scope: one platform zone, one side, one AHU.
# The station has two symmetric zones (one per side). Both are identical —
# model one, results apply to both. The two zones are NOT thermally coupled.
#
# Convention:
#   [SOBOL]      — parameter swept in Sobol sensitivity analysis. Range noted.
#   [ASSUMPTION] — engineering estimate. No hard reference. Justified below.
#   [STUDY]      — value requires a dedicated engineering study to compute
#                  rigorously. Assumption taken explicitly because that study
#                  is out of scope for a public-data model.
# =============================================================================

import numpy as np

# -----------------------------------------------------------------------------
# 1. GEOMETRY — one platform zone (one side)
# -----------------------------------------------------------------------------

PLATFORM_LENGTH_M    = 55.0
# m — usable platform length.
# Basis: GPE ligne 15 uses MP14 rolling stock, 3-car formation.
# 3 cars × ~18.5 m/car ≈ 55.5 m → 55 m used.
# Source: SGP rolling stock technical note (MP14), referenced in
# Dossier d'Enquête Publique ligne 15 Sud (2015):
# https://www.societedugrandparis.fr/gpe/ligne/ligne-15-sud

PLATFORM_WIDTH_M     =  4.0
# m — usable waiting zone width, one side.
# Derived from occupancy: 250 p on 220 m² = 1.14 p/m² (Fruin LOS C).
# 220 m² / 55 m = 4.0 m. Width is a consequence of the occupancy
# assumption, not an independent measurement.
# [ASSUMPTION] — no architectural drawing used.
# [SOBOL — range 3.5–5.0 m]

PLATFORM_HEIGHT_M    =  4.2
# m — floor-to-ceiling height of the platform zone.
# [ASSUMPTION] — no GPE architectural drawing used.
# Typical value for modern underground cut-and-cover metro stations,
# consistent with published cross-sections of Paris M14 extensions and
# London Jubilee Line deep-level stations.
# A rigorous value requires the station architectural plans — a civil
# engineering deliverable not available in public data for any GPE station.
# [SOBOL — range 3.8–5.0 m]

FACADE_HEIGHT_M      =  2.8
# m — PSD glass wall height (floor slab to top of glass panel).
# [ASSUMPTION] — no GPE drawing used.
# Constrained by: train door height (~2.1 m) + structural clearance above.
# Standard for modern metro PSD assemblies (Paris M14 reference).
# Rigorous value requires the PSD architectural/structural drawings —
# a manufacturer and civil engineering deliverable not publicly available.
# [SOBOL — range 2.5–3.2 m]

# Derived geometry
PLATFORM_AREA_M2     = PLATFORM_LENGTH_M * PLATFORM_WIDTH_M
# 220.0 m²

V_PLATFORM_M3        = PLATFORM_LENGTH_M * PLATFORM_WIDTH_M * PLATFORM_HEIGHT_M
# 924.0 m³

A_FACADE_M2          = PLATFORM_LENGTH_M * FACADE_HEIGHT_M
# 154.0 m² — PSD glass wall, platform-tunnel interface

A_SOIL_M2            = (
      PLATFORM_LENGTH_M * PLATFORM_HEIGHT_M          # back wall (long):  231.0
    + 2 * PLATFORM_WIDTH_M * PLATFORM_HEIGHT_M       # two end walls:      33.6
    + PLATFORM_LENGTH_M * PLATFORM_WIDTH_M           # ceiling:           220.0
    + PLATFORM_LENGTH_M * PLATFORM_WIDTH_M           # floor:             220.0
)
# 704.6 m² — concrete surfaces in contact with soil


# -----------------------------------------------------------------------------
# 2. THERMAL ENVELOPE — split into two boundary conditions
# -----------------------------------------------------------------------------
# The station is fully buried. T_ext does NOT appear directly as an envelope
# boundary condition. Two distinct thermal interfaces:
#   Facade (PSD glass wall) → tunnel air at T_tun   [high U, dynamic V_inf]
#   Structure (concrete box) → stable soil at T_soil [low U, quasi-steady]

# --- 2a. Facade: PSD glass wall + metal frame assembly ---

U_FACADE_W_M2K       =  7.0
# W/m²K — thermal transmittance of the PSD facade assembly.
#
# Derivation (traceable chain):
#   R_si  = 0.13 m²K/W  — internal still air, platform side
#             [ASHRAE Fundamentals 2017, Ch.26, Table 1]
#   R_glass = 0.012 m²K/W — 12 mm tempered/laminated glass, λ = 1.0 W/mK
#             [EN ISO 10456:2007, Table 1]
#   R_se  = 0.04 m²K/W  — tunnel side, moving air (piston effect)
#             [ASHRAE Fundamentals 2017, Ch.26, Table 1, forced convection]
#   → U_glass_only = 1 / (0.13 + 0.012 + 0.04) = 5.5 W/m²K
#
#   Metal frame effect (uninsulated aluminum/steel frames):
#   ASHRAE 90.1-2019, Table A2.3 + EN ISO 10077-1:2017:
#   "Uninsulated metal frame + single glazing" → assembly U = 6.0–8.0 W/m²K.
#   → U = 7.0 W/m²K taken as mid-range, consistent with heavy-duty
#     industrial PSD specification.
#
# A precise value requires a certified laboratory thermal test of the actual
# PSD product (EN ISO 12567 test protocol) — a manufacturer qualification
# deliverable, not publicly available for any GPE station PSD supplier.
# [SOBOL — range 5.5–8.5 W/m²K]

UA_FACADE_W_K        = U_FACADE_W_M2K * A_FACADE_M2
# 1078.0 W/K — static conduction term (doors closed)
# The total effective facade conductance is: UA_FACADE + RHO_CP_AIR * V_inf
# where V_inf is dynamic (function of train frequency). See §6 and occupancy.py.

# --- 2b. Structural envelope: concrete walls, ceiling, floor ---

U_SOIL_W_M2K         =  0.41
# W/m²K — effective thermal transmittance of buried concrete structure.
#
# Derivation:
#   R_si      = 0.13  m²K/W  [ASHRAE Fundamentals 2017, Ch.26]
#   R_concrete = 0.29 m²K/W  — 0.5 m wall, λ = 1.7 W/mK
#                               [EN ISO 10456:2007, reinforced concrete]
#   R_soil    = 2.00  m²K/W  — effective path d_eff = 3.0 m, k_soil = 1.5 W/mK
#                               (damp sandy clay, Paris basin, above water table)
#                               [ISO 13370:2017 §6; ASHRAE Fundamentals 2017, Ch.25]
#   → U_s = 1 / (0.13 + 0.29 + 2.00) = 1 / 2.42 = 0.41 W/m²K
#
# R_soil sensitivity to soil moisture:
#   Dry soil     (k ≈ 1.0): R_soil = 3.0 → U_s = 0.35 W/m²K
#   Saturated    (k ≈ 2.5): R_soil = 1.2 → U_s = 0.59 W/m²K
#
# A rigorous U_s requires a geotechnical site investigation: soil boring logs,
# in-situ thermal conductivity measurements (needle probe or TRT), and
# groundwater level monitoring over at least one annual cycle. This is a
# dedicated geotechnical engineering study — not available from public data
# for any specific GPE station site.
# [SOBOL — range 0.35–0.59 W/m²K]

UA_SOIL_W_K          = U_SOIL_W_M2K * A_SOIL_M2
# 289.0 W/K


# -----------------------------------------------------------------------------
# 3. BOUNDARY TEMPERATURES
# -----------------------------------------------------------------------------

T_SOIL_C             = 15.0
# °C — undisturbed ground temperature at depth > 5 m, Île-de-France.
# Measured range: 13–15°C. Upper bound used — conservative (maximises
# summer cooling load, worst-case design).
# Source: BRGM (Bureau de Recherches Géologiques et Minières),
# subsurface temperature atlas of metropolitan France.
# https://www.brgm.fr/fr/enjeu/geothermie
# [SOBOL — range 13–15°C]

T_TUN_OFFSET_C       =  5.0
# °C — tunnel air excess temperature above outdoor air: T_tun = T_ext + 5.
# Physical origin: heat dissipated into tunnel air by train braking
# (resistive or regenerative losses), traction motors, and passenger body
# heat inside train cars. Range from SES (Subway Environment Simulation)
# literature: 3–8°C depending on headway, fleet thermal profile,
# and tunnel ventilation effectiveness.
# [ASSUMPTION] — no public temperature measurement for GPE ligne 15
# tunnels (line not yet in commercial operation).
# [STUDY] — a rigorous T_tun requires a full SES aeraulic and thermodynamic
# simulation of the tunnel-station system, incorporating train speed profiles,
# motor heat rejection curves (from traction system specs), braking energy
# balance, and tunnel ventilation capacity. This is a dedicated engineering
# study (typically 6–12 months, conducted by the tunnel ventilation specialist
# during detailed design phase). Entirely out of scope for a public-data model.
# [SOBOL — range 3–8°C, HIGH PRIORITY: dominant uncertainty in T_tun]


# -----------------------------------------------------------------------------
# 4. THERMAL CAPACITANCE
# -----------------------------------------------------------------------------

RHO_AIR_KG_M3        =  1.2
CP_AIR_J_KG_K        = 1005.0
# Standard dry air at ~20°C, 1 atm.
# Source: ASHRAE Fundamentals 2017, Chapter 1.

RHO_CONC_KG_M3       = 2300.0
CP_CONC_J_KG_K       =  880.0
# Reinforced concrete.
# Source: EN 1991-1-1:2002 (density); EN ISO 10456:2007 (specific heat).

D_CONC_EFF_M         =  0.12
# m — effective concrete depth for active thermal storage.
# Physical basis: thermal penetration depth δ = sqrt(α × t),
#   α = λ/(ρ×cp) = 1.7 / (2300 × 880) ≈ 8.4×10⁻⁷ m²/s
#   Over 1h:  δ ≈ 0.055 m
#   Over 4h:  δ ≈ 0.11 m
#   Over 24h: δ ≈ 0.27 m
# → 10–15 cm is consistent with diurnal (1–4h) station temperature cycles.
# Mid-range 12 cm used.
# Reference convention: CIBSE Guide A (2015), Chapter 3 (thermal response
# and dynamic simulation).
# [ASSUMPTION] — not derived from measurement. The active depth depends on
# actual thermal history and surface boundary conditions.
# Rigorous determination requires either: (a) in-situ temperature sensors
# at multiple depths in the structural concrete, or (b) a 2D/3D finite-element
# heat conduction model of the station box. Both are dedicated structural
# thermal engineering studies out of scope here.
# [SOBOL — range 0.05–0.20 m]

C_AIR_J_K            = RHO_AIR_KG_M3  * CP_AIR_J_KG_K  * V_PLATFORM_M3
# 1.115e6 J/K

C_CONC_J_K           = RHO_CONC_KG_M3 * CP_CONC_J_KG_K * A_SOIL_M2 * D_CONC_EFF_M
# 1.718e8 J/K

C_TOTAL_J_K          = C_AIR_J_K + C_CONC_J_K
# 1.729e8 J/K
# Concrete dominates by 99%. Air mass is negligible for a buried concrete box.
# τ_static = C_TOTAL / (UA_FACADE + UA_SOIL) = 1.729e8 / 1367 ≈ 35 h
# The structural mass heats and cools very slowly.
# The AHU controls air temperature; the structure drifts slowly.


# -----------------------------------------------------------------------------
# 5. OCCUPANCY
# -----------------------------------------------------------------------------

PEOPLE_PEAK          = 250
# persons — peak occupancy, one platform side (one zone).
# Basis: Fruin LOS C, 220 m² platform area.
#   250 / 220 = 1.14 p/m² — LOS B/C boundary (busy, operationally acceptable).
# Source: Fruin, J.J. (1971). Designing for Pedestrians: A Level-of-Service
# Concept. Port Authority of New York.
# Referenced in: Metrolinx DS-12 Pedestrian Flow Modelling Design Standard,
# v1.0, May 2024, Section 5.2, Table 5-1.
# URL: https://assets.metrolinx.com/image/upload/v1738255845/Documents/
#      Engineering/DS-12_Pedestrian_Flow_Modelling_Design_Standard_v1.0.pdf
# [SOBOL — range 150–350 persons (LOS B to LOS D)]

WATTS_SENSIBLE_PP    =  75.0
# W/person — sensible heat, standing / slow walking at 24°C ambient.
# Source: ASHRAE Handbook of Fundamentals 2013, Chapter 18, Table 1.
# Row: "Standing, light work; walking" — Sensible = 75 W, Total = 160 W.
# [SOBOL — range 70–90 W]

WATTS_LATENT_PP      =  60.0
# W/person — latent heat. Humidity submodel only. NOT in thermal ODE.
# Source: same ASHRAE table. Range 55–70 W, mid value used.

BASELINE_W           = 5_000.0
# W — continuous platform load: LED lighting only.
# Escalators, gates, screens are concourse items — excluded.
# Derivation:
#   SEAM4US Barcelona (Ansuini et al. 2013), Table 1:
#   264 T8 lamps × 36 W = 9504 W fluorescent baseline.
#   Source: https://www.irbnet.de/daten/iconda/CIB_DC26414.pdf
#   LED correction: ×0.50 → 4752 W → 5000 W.
# [SOBOL — range 3000–8000 W]


# -----------------------------------------------------------------------------
# 6. INFILTRATION — PSD door openings
# -----------------------------------------------------------------------------
# Two mechanisms drive infiltration even with full-height PSDs:
#
# A) DOOR-OPEN BULK EXCHANGE (trains stopped, PSDs open)
#    All 12 PSDs per side open simultaneously for ~20 s per stop. The
#    piston effect of the arriving train has pressurised the tunnel ahead
#    of it. On door opening, the pressure differential drives a transient
#    air exchange — tunnel air enters the platform volume.
#
# B) GAP INFILTRATION (PSDs closed)
#    PSD assemblies have imperfect sealing: gaps at door edges, top/bottom
#    frame joints, and service penetrations. Train piston pressure waves
#    transiently push tunnel air through these gaps. Platform overpressure
#    (see §7) keeps steady-state flow direction outward, but does not
#    eliminate transient spikes during each train arrival.
#
# Method chosen: Exchange Efficiency
#   V_cycle = η × V_platform  (m³ exchanged per train stop)
#   V̇_inf   = V_cycle / T_headway  (time-averaged, m³/s)
#
# Why NOT the velocity method (A_open × v_air × Δt)?
#   The velocity method requires v_air: the air velocity through the open
#   door area during the piston pressure event. This cannot be estimated
#   without a full aeraulic simulation of the tunnel-station interface —
#   including train speed profile, tunnel cross-section geometry, PSD
#   leakage coefficient, and platform overpressure level. This is a
#   dedicated SES (Subway Environment Simulation) engineering study,
#   typically conducted by the tunnel ventilation specialist over several
#   months of detailed design. It is explicitly out of scope for a
#   public-data model. The exchange efficiency method sidesteps this by
#   using η, which is better constrained by literature measurements.
#
# η = 0.15 — fraction of platform air exchanged per train stop.
# Range from PSD subway measurement literature: 0.10–0.25.
# Search: "Experimental study on air exchange coefficient of subway
#          stations with platform screen doors" (ScienceDirect).
# [ASSUMPTION] — no GPE-specific or Paris Metro measurement found.
# [STUDY] — rigorous η requires either:
#   (a) tracer gas measurements during live train operations, or
#   (b) full CFD/SES simulation of the tunnel-station system.
#   Both are dedicated engineering studies out of scope here.
# [SOBOL — range 0.10–0.25, HIGH PRIORITY]

ETA_INF              = 0.15
V_CYCLE_M3           = ETA_INF * V_PLATFORM_M3
# 138.6 m³ exchanged per train stop

# --- Train headways — GPE ligne 15 Ouest ---
# Source (peak): Dossier d'Enquête d'Utilité Publique (DUP), ligne 15,
# Pièce G — Caractéristiques de l'infrastructure, Section 3 (performances).
# DUP publicly available: https://www.societedugrandparis.fr/gpe/ligne/ligne-15-sud
# Committed peak headway: 2 min (base model).
# Ultimate capacity: 1.5 min (SGP capacity planning communications, 2022) —
# used as Sobol upper bound only.
#
# Off-peak / night headways: [ASSUMPTION] — DUP gives peak headway only.
# Derived from typical metro scheduling practice.
# Night: ligne 15 has no all-night service in current planning.

HEADWAY_PEAK_S       = 120.0    # s — 2 min, 30 trains/h  [DUP base]
HEADWAY_OFFPEAK_S    = 240.0    # s — 4 min, 15 trains/h  [ASSUMPTION]
HEADWAY_NIGHT_S      = np.inf   # no service

HEADWAY_PEAK_SOBOL_MIN_S = 90.0  # s — 1.5 min, 40 trains/h [Sobol upper bound]

# V̇_inf is computed dynamically — see occupancy.v_inf_m3s(hour, day_type).


# -----------------------------------------------------------------------------
# 7. OVERPRESSURE — minimum platform airflow
# -----------------------------------------------------------------------------
# Rationale for maintaining overpressure even with full-height PSDs:
#
#   1. GAP INFILTRATION BARRIER: Imperfect PSD sealing means every gap is a
#      potential infiltration point. Platform overpressure ensures net air
#      flow direction is always platform→tunnel through all gaps, preventing
#      hot, humid, dusty tunnel air from entering through edges and joints.
#
#   2. GLASS STRUCTURAL STRESS: A pressure differential across the PSD facade
#      creates mechanical loading on glass panels and frame connections.
#      Overpressure keeps the differential within the structural design limit
#      (platform slightly above tunnel, not reverse — reverse loading is the
#      design failure mode).
#
#   3. SMOKE CONTROL (fire/emergency): Regulations require that in case of
#      a fire event in the tunnel, smoke must be directed away from passengers.
#      Platform overpressure is part of the smoke control strategy — it pushes
#      combustion products toward tunnel extraction points.
#
# Sizing:
# [ASSUMPTION] [STUDY] — the minimum airflow needed to maintain a target
# overpressure (typically 10–30 Pa above tunnel pressure) depends on:
# total PSD leakage area (sum of all gap cross-sections), platform volume,
# and AHU pressure capability. This requires a full aeraulic pressure balance
# study of the platform-PSD-tunnel system — a dedicated ventilation engineering
# deliverable (part of the SES study package and fire safety engineering report
# for the station). This study is conducted by specialist ventilation engineers
# over the detailed design phase and is not available from public data.
# Value of 2500 m³/h per zone is an engineering estimate consistent with
# typical practice for PSD-equipped underground metro stations.
# [SOBOL — range 1500–4000 m³/h]

AIRFLOW_OVERPRESSURE_M3H = 2_500.0   # m³/h — per zone (one side, one AHU)


# -----------------------------------------------------------------------------
# 8. AIRFLOW SIZING — one zone (one AHU, one platform side)
# -----------------------------------------------------------------------------

AIRFLOW_PER_PERSON_M3H = 25.0
# m³/h/person — minimum fresh air per occupant.
# Source: EN 16798-1:2019, Category II (normal expectation), high-density
# transit occupancy. Consistent with ERP type GA (French code, gares).

# Additive logic — Q_total = Q_overpressure + Q_occupancy:
# The two requirements are physically independent:
#   Overpressure demand: present at all times regardless of headcount.
#   Hygiene demand: scales with occupancy, independent of pressure.
# They serve different physical purposes and cannot substitute for each other.
# Taking max() would be wrong: at zero occupancy, Q_occupancy = 0 but
# Q_overpressure must still be delivered in full.

AIRFLOW_MIN_M3H = AIRFLOW_OVERPRESSURE_M3H
# 2500 m³/h — zero occupancy, overpressure maintained continuously.

AIRFLOW_MAX_M3H = (AIRFLOW_OVERPRESSURE_M3H
                   + PEOPLE_PEAK * AIRFLOW_PER_PERSON_M3H) * 1.10
# (2500 + 250×25) × 1.10 = 8750 × 1.10 = 9625 m³/h
# 10% safety margin on peak demand.


# -----------------------------------------------------------------------------
# 9. REGULATION SETPOINT BOUNDARIES
# Documented fully in regulation.py. Values reproduced here for reference.
# Do NOT duplicate control logic — regulation.py is the authoritative source.
# -----------------------------------------------------------------------------

T_ANTIFREEZE_C       =  5.0   # °C — anti-freeze protection floor (T_in minimum)
T_HEAT_LINEAR_LOW_C  = -1.0   # °C — T_ext threshold below which anti-freeze mode
T_HEAT_LINEAR_HIGH_C =  6.0   # °C — T_ext threshold above which fixed 12°C target
T_HEAT_FIXED_C       = 12.0   # °C — fixed heating target (T_ext 6–12°C zone)
T_DEAD_LOW_C         = 12.0   # °C — dead band lower boundary (T_ext)
T_DEAD_HIGH_C        = 26.0   # °C — dead band upper boundary (T_ext)
T_COOL_FIXED_C       = 26.0   # °C — fixed cooling target (T_ext 26–31°C zone)
T_COOL_LINEAR_HIGH_C = 31.0   # °C — T_ext above which linear cooling offset
T_BLOW_COOL_C   = 15.0   # °C — fixed supply air temp, cooling mode [ASSUMPTION]
T_BLOW_HEAT_C   = 30.0   # °C — fixed supply air temp, heating mode [ASSUMPTION]
# Cooling: 13–16°C typical (below avoids condensation). Heating: 28–35°C typical.
# [SOBOL — cool range 13–17, heat range 28–35]

# -----------------------------------------------------------------------------
# 10. AIR PROPERTIES — convenience aliases
# -----------------------------------------------------------------------------

CP_AIR_KJ_KG_K       = CP_AIR_J_KG_K / 1000.0
# kJ/(kg·K) — used in HVAC power calculations

RHO_CP_AIR_J_M3_K   = RHO_AIR_KG_M3 * CP_AIR_J_KG_K
# 1206 J/(m³·K) — volumetric heat capacity of air
# Used in infiltration heat flux: Q_inf = RHO_CP_AIR × V̇_inf × (T_tun - T_in)
# and in HVAC power:              Q_hvac = RHO_CP_AIR × Q_air_m3s × ΔT
