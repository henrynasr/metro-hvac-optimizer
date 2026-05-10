# =============================================================================
# constants.py — Energy Twin: single source of truth for all model parameters.
# One platform zone, one side, one AHU. See parameters.md for full derivations.
# Convention: [SOBOL min–max] marks parameters swept in sensitivity analysis.
# =============================================================================

import numpy as np

# -----------------------------------------------------------------------------
# 1. GEOMETRY
# -----------------------------------------------------------------------------

PLATFORM_LENGTH_M    = 55.0    # m  — MP14 3-car train length
PLATFORM_WIDTH_M     =  4.0    # m  — derived from 250p at Fruin LOS C. [ASSUMPTION] [SOBOL 3.5–5.0]
PLATFORM_HEIGHT_M    =  4.2    # m  — typical cut-and-cover metro. [ASSUMPTION] [SOBOL 3.8–5.0]
FACADE_HEIGHT_M      =  2.8    # m  — PSD glass wall height. [ASSUMPTION] [SOBOL 2.5–3.2]

PLATFORM_AREA_M2     = PLATFORM_LENGTH_M * PLATFORM_WIDTH_M          # 220.0 m²
V_PLATFORM_M3        = PLATFORM_LENGTH_M * PLATFORM_WIDTH_M * PLATFORM_HEIGHT_M  # 924.0 m³
A_FACADE_M2          = PLATFORM_LENGTH_M * FACADE_HEIGHT_M            # 154.0 m² — PSD wall
A_STAIR_WALL_M2      = PLATFORM_WIDTH_M * PLATFORM_HEIGHT_M            # 16.8 m² — stair entry wall (open, excluded from soil)
A_SOIL_M2            = (
      PLATFORM_LENGTH_M * PLATFORM_HEIGHT_M        # back wall:   231.0
    + 1 * PLATFORM_WIDTH_M * PLATFORM_HEIGHT_M     # one end wall:  16.8 (stair side excluded)
    + PLATFORM_LENGTH_M * PLATFORM_WIDTH_M         # ceiling:     220.0
    + PLATFORM_LENGTH_M * PLATFORM_WIDTH_M         # floor:       220.0
)                                                                      # 687.8 m²

# -----------------------------------------------------------------------------
# 2. THERMAL ENVELOPE
# -----------------------------------------------------------------------------

U_FACADE_W_M2K       =  7.0    # W/m²K — PSD glass + metal frame assembly. [SOBOL 5.5–8.5]
UA_FACADE_W_K        = U_FACADE_W_M2K * A_FACADE_M2                   # 1078.0 W/K (static, doors closed)

A_TUN_WALL_M2        = (PLATFORM_HEIGHT_M - FACADE_HEIGHT_M) * PLATFORM_LENGTH_M  # 77.0 m² — concrete wall above PSD, faces tunnel
U_TUN_WALL_W_M2K     =  2.5    # W/m²K — concrete wall, tunnel side. R_si=0.13 + R_conc(0.3m)=0.18 + R_se(forced)=0.04 -> U=2.9, take 2.5 conservative. [ASSUMPTION] [SOBOL 2.0–3.5]
UA_TUN_WALL_W_K      = U_TUN_WALL_W_M2K * A_TUN_WALL_M2                  # 192.5 W/K — same boundary as facade (T_tun)

U_SOIL_W_M2K         =  0.41   # W/m²K — buried concrete box (R_si + R_concrete + R_soil). [SOBOL 0.35–0.59]
UA_SOIL_W_K          = U_SOIL_W_M2K * A_SOIL_M2                       # 289.0 W/K

# -----------------------------------------------------------------------------
# 3. BOUNDARY TEMPERATURES
# -----------------------------------------------------------------------------

T_SOIL_C             = 15.0    # °C — stable ground temp, Île-de-France (BRGM). [SOBOL 13–17]
T_TUN_OFFSET_DAY_C   = 10.0    # °C — T_tun = T_ext + 10 during service (05h–23h). [SOBOL 5–15]
T_TUN_OFFSET_NIGHT_C =  5.0    # °C — T_tun = T_ext + 5 at night (23h–05h).       [SOBOL 3–8]

# -----------------------------------------------------------------------------
# 4. THERMAL CAPACITANCE
# -----------------------------------------------------------------------------

RHO_AIR_KG_M3        =  1.2    # kg/m³
CP_AIR_J_KG_K        = 1005.0  # J/(kg·K)
RHO_CONC_KG_M3       = 2300.0  # kg/m³
CP_CONC_J_KG_K       =  880.0  # J/(kg·K)

D_CONC_EFF_M         =  0.12   # m — effective concrete depth for thermal storage (δ=√αt, ~4h cycle). [SOBOL 0.05–0.20]

C_AIR_J_K            = RHO_AIR_KG_M3  * CP_AIR_J_KG_K  * V_PLATFORM_M3           # 1.115e6 J/K
C_CONC_J_K           = RHO_CONC_KG_M3 * CP_CONC_J_KG_K * A_SOIL_M2 * D_CONC_EFF_M  # 1.718e8 J/K
C_TOTAL_J_K          = C_AIR_J_K + C_CONC_J_K                                     # 1.729e8 J/K  τ≈35h

CP_GLYCOL_J_KG_K     = 3800.0  # J/(kg·K) — 30% propylene glycol. [ASSUMPTION]
RHO_GLYCOL_KG_M3     = 1045.0  # kg/m³    — 30% propylene glycol. [ASSUMPTION]

# -----------------------------------------------------------------------------
# 5. OCCUPANCY
# -----------------------------------------------------------------------------

PEOPLE_PEAK          = 250     # persons — Fruin LOS C on 220 m². [SOBOL 150–350]
WATTS_SENSIBLE_PP    =  75.0   # W/person — sensible heat, standing/slow walk (ASHRAE 2013 Ch.18). [SOBOL 70–90]
WATTS_LATENT_PP      =  60.0   # W/person — latent heat (psychrometric layer only, not in ODE).
BASELINE_W           = 500.0   # W — LED strips + arrival screens, corridor only. [SOBOL 300–800]

# -----------------------------------------------------------------------------
# 6. INFILTRATION — PSD door openings (exchange efficiency method)
# -----------------------------------------------------------------------------

ETA_INF              = 0.15    # — fraction of platform air exchanged per train stop. [SOBOL 0.10–0.25]
V_CYCLE_M3           = ETA_INF * V_PLATFORM_M3                        # 138.6 m³/stop

HEADWAY_PEAK_S       = 120.0   # s — 2 min peak (GPE DUP).
HEADWAY_OFFPEAK_S    = 240.0   # s — 4 min off-peak. [ASSUMPTION]
HEADWAY_NIGHT_S      = np.inf  # s — no service.
HEADWAY_PEAK_SOBOL_MIN_S = 90.0  # s — 1.5 min upper bound for Sobol.

# -----------------------------------------------------------------------------
# 6b. STAIRCASE FRESH AIR — passive outdoor air via open stair entrance
# -----------------------------------------------------------------------------
# Q_stair = V_AIR_STAIR * A_STAIR  (fixed flow, constant)
# Heat term: RHO_CP_AIR * Q_STAIR_M3S * (T_ext - T_in)

A_STAIR_M2           =  2.5   # m² — stair free area. [ASSUMPTION] [SOBOL 2.0–3.5]
V_AIR_STAIR_MS       =  0.5   # m/s — mean air velocity in staircase. [SOBOL 0.3–0.7]
Q_STAIR_M3S          = V_AIR_STAIR_MS * A_STAIR_M2                    # 1.25 m³/s = 4500 m³/h

# -----------------------------------------------------------------------------
# 7. OVERPRESSURE — minimum AHU airflow
# -----------------------------------------------------------------------------

AIRFLOW_OVERPRESSURE_M3H = 2_500.0  # m³/h — gap barrier + smoke control. [ASSUMPTION] [SOBOL 1500–4000]

# -----------------------------------------------------------------------------
# 8. AIRFLOW SIZING
# -----------------------------------------------------------------------------

AIRFLOW_PER_PERSON_M3H = 25.0   # m³/h/person — EN 16798-1:2019 Cat. II minimum.
AIRFLOW_MIN_M3H = AIRFLOW_OVERPRESSURE_M3H                            # 2500 m³/h
AIRFLOW_MAX_M3H = (AIRFLOW_OVERPRESSURE_M3H + PEOPLE_PEAK * AIRFLOW_PER_PERSON_M3H) * 1.10  # 9625 m³/h

# -----------------------------------------------------------------------------
# 9. REGULATION SETPOINTS
# -----------------------------------------------------------------------------

T_EXT_ANTIFREEZE_C   =  0.0    # °C — anti-freeze mode on
T_ANTIFREEZE_C       =  5.0    # °C — anti-freeze floor
T_HEAT_FIXED_C       = 15.0    # °C — heating target (T_ext < T_DEAD_LOW). Source: metro platform comfort studies.
T_DEAD_LOW_C         = 15.0    # °C — dead band lower boundary (T_ext)
T_DEAD_HIGH_C        = 26.0    # °C — dead band upper boundary (T_ext)
T_COOL_FIXED_C       = 26.0    # °C — cooling target (T_ext > T_DEAD_HIGH). Source: metro platform comfort studies.
T_COOL_DELTA_C       =  6.0    # °C — cooling target = T_ext - 6 (e.g. 34°C when T_ext=40). [SOBOL 4–8]
T_BLOW_COOL_C        = 15.0    # °C — AHU supply, cooling mode. [ASSUMPTION] [SOBOL 13–17]
T_BLOW_HEAT_C        = 30.0    # °C — AHU supply, heating mode. [ASSUMPTION] [SOBOL 28–35]

# -----------------------------------------------------------------------------
# 10. AIR PROPERTIES
# -----------------------------------------------------------------------------

CP_AIR_KJ_KG_K       = CP_AIR_J_KG_K / 1000.0   # kJ/(kg·K)
RHO_CP_AIR_J_M3_K    = RHO_AIR_KG_M3 * CP_AIR_J_KG_K  # 1206 J/(m³·K)

# -----------------------------------------------------------------------------
# 11. WATER REGIME
# -----------------------------------------------------------------------------

T_HW_EXT_LOW_C       = -7.0    # °C — T_ext anchor, hot water max supply
T_HW_EXT_HIGH_C      = 15.0    # °C — hot water shutoff
T_HW_EXT_HYST_C      = 13.0    # °C — hot water restart threshold
T_HW_SUPPLY_MAX      = 50.0    # °C — hot water supply at T_ext = -7
T_HW_SUPPLY_MIN      = 35.0    # °C — hot water supply at T_ext = 12
T_HW_RETURN_MAX      = 45.0    # °C
T_HW_RETURN_MIN      = 30.0    # °C

T_CW_EXT_LOW_C       = 26.0    # °C — cold water shutoff
T_CW_EXT_HIGH_C      = 31.0    # °C — T_ext anchor, cold water min supply
T_CW_EXT_HYST_C      = 27.0    # °C — cold water restart threshold
T_CW_SUPPLY_MIN      =  8.0    # °C — cold water supply at T_ext = 31
T_CW_SUPPLY_MAX      = 12.0    # °C — cold water supply at T_ext = 26
T_CW_RETURN_MIN      = 14.0    # °C
T_CW_RETURN_MAX      = 18.0    # °C

DT_WATER_HEAT_K      = T_HW_SUPPLY_MAX - T_HW_RETURN_MAX   # 5 K — heating circuit ΔT
DT_WATER_COOL_K      = T_CW_RETURN_MIN - T_CW_SUPPLY_MIN   # 6 K — cooling circuit ΔT

# -----------------------------------------------------------------------------
# 12. AHU AIR MIX
# -----------------------------------------------------------------------------

FRAC_RETURN_AIR      = 0.70    # — 70% return / 30% outdoor. Industry practice (Seoul SCAP, Delhi Metro).

# -----------------------------------------------------------------------------
# 13. HUMIDITY TARGETS
# -----------------------------------------------------------------------------

RH_TARGET_LOW        = 0.40    # — comfort lower bound (ASHRAE 55)
RH_TARGET_HIGH       = 0.60    # — comfort upper bound / dehumidification target
RH_ALERT_LOW         = 0.30    # — dryness alert
RH_ALERT_HIGH        = 0.70    # — condensation risk threshold
RH_STRUCTURAL        = 0.80    # — ASHRAE 160 mold criterion, 30-day avg

# -----------------------------------------------------------------------------
# 14. Emissions 
# -----------------------------------------------------------------------------

COP_COOL  = 5.5   # water-cooled chiller             [SOBOL 4.0–6.5]
COP_HEAT  = 4.0   # water-water heat pump, 50°C supply [SOBOL 3.0–5.5]
ETA_FAN   = 0.60  # fan + motor + VFD combined        [SOBOL 0.40–0.70]
DP_AHU_PA = 800   # total static pressure, AHU+ducts  [SOBOL 500–1200]
ELEC_PRICE_EUR_KWH = 0.17   # €/kWh — RATP industrial tariff [ASSUMPTION]