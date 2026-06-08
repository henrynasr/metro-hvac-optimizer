"""
Fault Detection & Diagnostics — T_ext sensor drift detector.

Idea:
  - "Reality": run the model on CLEAN weather -> gives the true T_in
    (what the wall probe would read).
  - "Sensor": same weather but +2 C after July 1 -> the corrupted
    T_ext stream the detector reads, one hour at a time.
  - Each hour: solve the ODE forward ONE step on the corrupted T_ext,
    compare predicted T_in to the true T_in. If T_ext is honest they
    match; once it drifts, a gap opens.
  - Two-sided CUSUM on that gap. Halt the instant it trips.
"""

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from utils import load_data
from occupancy import load_profiles, build_Q_array
from regulation import dT_dt

# --------------------------------------------------------------------
# 1. Inputs (same as simulation.py — one full year 2024)
# --------------------------------------------------------------------
SECONDS_PER_HOUR = 3600
FAULT_START = "2024-07-01"   # sensor starts drifting here
FAULT_BIAS = 2.0             # +2 C bias
K = 0.05                      # CUSUM slack: ignore noise below 0.5 C
H = 3.0                      # CUSUM alarm threshold (calibrate on clean run)

df = load_data("data/raw/paris_weather.csv")          # adjust path if needed
df = df.loc["2024-01-01":"2024-12-31"]
dates = df.index

t_array = np.arange(len(df)) * SECONDS_PER_HOUR
T_ext_clean = df["temperature_2m"].values

profiles = load_profiles("data/raw/Defense_Occupation_Normalised.xlsx")
Q_int_array, n_people_array = build_Q_array(dates, profiles)

# --------------------------------------------------------------------
# 2. Build the two T_ext streams
# --------------------------------------------------------------------
# Corrupted sensor: +2 C after July 1
T_ext_sensor = T_ext_clean.copy()
fault_mask = dates >= FAULT_START
T_ext_sensor[fault_mask] += FAULT_BIAS

# --------------------------------------------------------------------
# 3. Helper: solve the ODE for ONE hour, return T_in at end of hour
# --------------------------------------------------------------------

water_state = {"heating": False, "cooling": False}

def solve_one_hour(T_in_start, hour_index, T_ext_stream):
    n = len(T_ext_stream)
    if hour_index + 2 <= n:
        t_ext_slice    = T_ext_stream[hour_index:hour_index+2]
        q_int_slice    = Q_int_array[hour_index:hour_index+2]
        n_ppl_slice    = n_people_array[hour_index:hour_index+2]
        dates_slice    = dates[hour_index:hour_index+2]
    else:
        t_ext_slice    = T_ext_stream[hour_index:hour_index+1].repeat(2)
        q_int_slice    = Q_int_array[hour_index:hour_index+1].repeat(2)
        n_ppl_slice    = n_people_array[hour_index:hour_index+1].repeat(2)
        dates_slice    = dates[hour_index:hour_index+1].repeat(2)

    sol = solve_ivp(
        dT_dt,
        t_span=(0, SECONDS_PER_HOUR),
        y0=[T_in_start],
        args=(
            np.array([0, SECONDS_PER_HOUR]),
            t_ext_slice,
            q_int_slice,
            n_ppl_slice,
            dates_slice,
            water_state,
        ),
        t_eval=[SECONDS_PER_HOUR],
        method="RK45",
        max_step=3600.0,
    )
    return sol.y[0][-1]

# --------------------------------------------------------------------
# 4. "Reality": true T_in from CLEAN weather (the wall probe)
# --------------------------------------------------------------------
T_in_measured = np.zeros(len(df))
T_in = T_ext_clean[0]                       # seed first hour
for t in range(len(df)):
    T_in_measured[t] = solve_one_hour(T_in, t, T_ext_clean)
    T_in = T_in_measured[t]

# --------------------------------------------------------------------
# 5. Streaming detector: read corrupted sensor hour by hour
# --------------------------------------------------------------------
S_hi = 0.0
S_lo = 0.0
T_in = T_in_measured[0]                      # start from real probe
alarm_hour = None

for t in range(len(df)):
    # predict T_in this hour using the SUSPECT T_ext
    T_in_pred = solve_one_hour(T_in, t, T_ext_sensor)

    # residual: what the probe really reads vs what the model expects
    resid = T_in_measured[t] - T_in_pred

    # two-sided CUSUM
    S_hi = max(0.0, S_hi + resid - K)
    S_lo = min(0.0, S_lo + resid + K)

    if t < 4344:
        print(t, round(resid, 3), round(S_hi, 2), round(S_lo, 2))

    if S_hi > H or S_lo < -H:
        alarm_hour = t
        break

    # next hour starts from the real probe reading
    T_in = T_in_measured[t]

# --------------------------------------------------------------------
# 6. Report
# --------------------------------------------------------------------
print("=" * 50)
print("  FDD — T_ext sensor drift detector")
print("=" * 50)
if alarm_hour is None:
    print("  No fault detected over the year.")
else:
    t_alarm = dates[alarm_hour]
    t_fault = pd.Timestamp(FAULT_START)
    lag_hours = (t_alarm - t_fault).total_seconds() / SECONDS_PER_HOUR
    lag_days = lag_hours / 24

    direction = "TOO HIGH" if S_hi > H else "TOO LOW"
    print(f"  ALARM fired")
    print(f"  When            : {t_alarm:%A %d %B %Y, %Hh}")
    print(f"  Fault injected  : {t_fault:%d %B %Y}")
    print(f"  Detection lag   : {lag_days:.1f} days ({lag_hours:.0f} h)")
    print(f"  Sensor reads    : {direction}")
print("=" * 50)
