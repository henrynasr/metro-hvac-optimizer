# =============================================================================
# constants.py — Energy Twin model parameters
# =============================================================================
# Single source of truth for all physical constants and model parameters.
# Import anywhere with: from constants import *
#
# All values validated Session 10 (2026-05-08). Sources documented in
# docs/parameters.md. Sobol analysis will quantify sensitivity to each.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. ZONE GEOMETRY — platform zone, one side, one AHU
# -----------------------------------------------------------------------------
PLATFORM_LENGTH_M   = 55.0   # m — platform usable length
PLATFORM_WIDTH_M    =  4.0   # m — usable waiting zone width (one side)
PLATFORM_HEIGHT_M   =  7.0   # m — floor to ceiling
PLATFORM_AREA_M2    = PLATFORM_LENGTH_M * PLATFORM_WIDTH_M          # 220 m²
AIR_VOLUME_M3       = PLATFORM_LENGTH_M * PLATFORM_WIDTH_M * PLATFORM_HEIGHT_M  # 1540 m³
ENVELOPE_AREA_M2    = (
    2 * (PLATFORM_LENGTH_M * PLATFORM_HEIGHT_M)   # two long walls
    + 2 * (PLATFORM_WIDTH_M  * PLATFORM_HEIGHT_M) # two end walls
    + 2 * (PLATFORM_LENGTH_M * PLATFORM_WIDTH_M)  # floor + ceiling
)  # 1266 m²

# -----------------------------------------------------------------------------
# 2. ENVELOPE — thermal conductance
# -----------------------------------------------------------------------------
# U = 0.5 W/m²K — engineering estimate, no direct source found for buried
# metro station walls. See docs/parameters.md §2 for full discussion.
U_WALL_W_M2K        =  0.5   # W/m²K — envelope U-value (estimate)
UA_W_K              = U_WALL_W_M2K * ENVELOPE_AREA_M2   # ~633 W/K → ~600 W/K

# -----------------------------------------------------------------------------
# 3. THERMAL CAPACITANCE — air mass only
# -----------------------------------------------------------------------------
# T_in represents air temperature only. Concrete mass excluded.
# Source: standard dry air properties at 20°C, 1 atm.
RHO_AIR_KG_M3       =  1.2    # kg/m³ — air density
CP_AIR_J_KG_K       = 1005.0  # J/(kg·K) — air specific heat
C_J_K               = RHO_AIR_KG_M3 * CP_AIR_J_KG_K * AIR_VOLUME_M3  # ~1.86e6 J/K

# -----------------------------------------------------------------------------
# 4. OCCUPANCY
# -----------------------------------------------------------------------------
# Peak: 250 persons/side → 1.14 p/m² → Fruin LOS C.
# Source: Fruin (1971) via Metrolinx DS-12 v1.0 (2024), Section 5.2.
PEOPLE_PEAK         = 250     # persons — peak occupancy, one platform side

# Sensible heat per occupant.
# Source: ASHRAE Fundamentals 2013, Chapter 18, Table 1.
# Activity: standing / light walking at 24°C room temperature.
WATTS_PER_PERSON    =  75.0   # W/person — sensible heat only (goes into ODE)
LATENT_W_PER_PERSON =  60.0   # W/person — latent heat (humidity model only, not ODE)

# -----------------------------------------------------------------------------
# 5. BASELINE EQUIPMENT LOAD — platform zone
# -----------------------------------------------------------------------------
# LED platform lighting only. Escalators/gates/screens are concourse items.
# Derived from SEAM4US Barcelona (Ansuini et al. 2013): 9.5 kW T8 fluorescent
# × 0.50 LED efficiency factor → ~4.8 kW → rounded to 5 kW.
BASELINE_W          = 5_000.0  # W — continuous equipment load (lighting)

# -----------------------------------------------------------------------------
# 6. AIR PROPERTIES (ventilation / HVAC power calculations)
# -----------------------------------------------------------------------------
CP_AIR_KJ_KG_K      = CP_AIR_J_KG_K / 1000   # kJ/(kg·K) — convenience alias

# -----------------------------------------------------------------------------
# 7. REGULATION — setpoint law boundaries
# Documented in regulation.py and README. Reproduced here for reference only.
# Do NOT duplicate logic — these are boundary values used in setpoint law.
# -----------------------------------------------------------------------------
T_ANTIFREEZE_C      =  5.0   # °C — minimum T_in (anti-freeze protection)
T_HEAT_FIXED_C      = 12.0   # °C — fixed heating target (T_ext 6–12°C zone)
T_DEAD_LOW_C        = 12.0   # °C — dead band lower boundary (T_ext)
T_DEAD_HIGH_C       = 26.0   # °C — dead band upper boundary (T_ext)
T_COOL_FIXED_C      = 26.0   # °C — fixed cooling target (T_ext 26–31°C zone)

# -----------------------------------------------------------------------------
# 8. AIRFLOW — sizing
# -----------------------------------------------------------------------------
# Sized on 25 m³/h/person, 2 sides, 2 AHUs. See README regulation section.
AIRFLOW_MIN_M3H     = 10_000.0  # m³/h — minimum (hygiene + overpressure)
AIRFLOW_MAX_M3H     = 22_000.0  # m³/h — peak (800 persons total, 2 AHUs)
AIRFLOW_PER_PERSON  =    25.0   # m³/h/person — regulatory minimum (EN 16798 / ERP GA)
