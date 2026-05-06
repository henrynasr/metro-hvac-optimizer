import numpy as np
import matplotlib.pyplot as plt


# System designed for these extreme temperatures 
T_ext_min = -7
T_ext_max = 32


# System temperature limits
T_anti_freeze = 5
delta_T_max_cooling = 5

def dT_dt(t, T, t_array, T_ext_array, Q_array, N_people_array, UA, C, ro=1.2, Cp=1000.0):
    """
    Slope of indoor temperature at instant t.

    Args:
        t (float): current time in seconds (since simulation start).
        T (float): current indoor temperature in °C.
        t_array (np.ndarray): time grid in seconds, shape (N,).
        T_ext_array (np.ndarray): outdoor temperature in °C at each
            point of t_array, shape (N,).
        Q_array (np.ndarray): internal heat load in W at each point
            of t_array, shape (N,).
        N_people_array (np.ndarray): Number of expected people at each
            point of t_array, shape (N,).
        UA (float): overall heat transfer coefficient * area, in W/K.
        C (float): lumped thermal capacitance, in J/K.
        ro (float): air density, in kg/m3
        Cp (flaot): air thermal capacity, in J/Kg.K

    Returns:
        float: dT/dt in °C/s.
    """
    T = float(T[0])
    T_ext = np.interp(t, t_array, T_ext_array)
    Q_internal = np.interp(t, t_array, Q_array)
    N_people = np.interp(t, t_array, N_people_array)

    Q_air = airflow_total(N_people)
    T_set = T_setpoint(T_ext)
    if np.isnan(T_set):
        Q_hvac = ro * Q_air * Cp * (T - T_ext) / 3600
    else: 
        Q_hvac = ro * Q_air * Cp * np.clip(T - T_set, -12, 5) / 3600 
    # air flow in hours not seconds
    # 12 = (+5-(-7)) worst case considered (anti-freezing)
    
    return (UA * (T_ext - T) + Q_internal-Q_hvac) / C

def build_Q_hvac_array(T, T_ext_array, N_people_array, ro = 1.2, Cp = 1000.0):
    n = len(T)
    Q_hvac_array = np.zeros(n)
    Q_vent_array = np.zeros(n)
    Q_heat_array = np.zeros(n)
    Q_cool_array = np.zeros(n)

    for i, T_in in enumerate(T):
        Q_air = airflow_total(N_people_array[i])
        T_ext = T_ext_array[i]
        T_set = T_setpoint(T_ext)

        if np.isnan(T_set):
            q = ro * Q_air * Cp * (T_in - T_ext_array[i]) / 3600
            Q_vent_array[i] = q
        else: 
            q = ro * Q_air * Cp * np.clip(T_in - T_set, T_ext_min - T_anti_freeze, delta_T_max_cooling) / 3600 
            if q<0:
                Q_heat_array[i] = q
            else:
                Q_cool_array[i] = q
        # air flow in hours not seconds
        # maximum delta T in heating = (T_ext_min - T_anti_freeze) worst case considered
        Q_hvac_array[i] = Q_vent_array[i] + Q_heat_array[i] + Q_cool_array[i]
    
    return Q_hvac_array, Q_heat_array, Q_cool_array, Q_vent_array

def T_setpoint(T_ext):
    """Return target T_in given T_ext. Returns None in dead band."""
    if T_ext < -1:
         return 5 # anti-freeze
    elif -1<= T_ext < 6:
        return T_ext + 6
    elif 6 <= T_ext < 12:
        return 12
    elif 12 <= T_ext <= 26:
        return np.nan # no control (dead band)
    elif 26 < T_ext <= 31:
        return 26
    else:
        return T_ext - 5

def airflow_total(n_people, n_peak=400, q_min=5000, q_max=11000):
    """Return total station airflow in m³/h."""
    return q_min + (q_max - q_min) * min(n_people / n_peak, 1.0)