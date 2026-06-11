import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go


# =============================================================================
# APP CONFIG — shared across all pages
# =============================================================================

st.set_page_config(page_title="Metro HVAC Optimizer", layout="wide")
st.title("Metro HVAC Optimizer")

page = st.sidebar.radio(
    "Page",
    ["Live Simulation", "Sobol Sensitivity", "Pareto Explorer"],
)


# =============================================================================
# PAGE 1 — LIVE SIMULATION
# Runs a full-year simulation with user-chosen parameters.
# =============================================================================

if page == "Live Simulation":

    st.caption("Runs a full-year simulation with your chosen parameters.")

    import constants
    import regulation
    import emissions as em_mod
    import humidity as hum_mod
    from simulation import run_simulation

    # --- Monkey-patch helpers (same pattern as pareto.py) ---
    def _patch(name, value):
        setattr(constants, name, value)
        for mod in (regulation, em_mod, hum_mod):
            if hasattr(mod, name):
                setattr(mod, name, value)

    def _patch_derived():
        ov = constants.AIRFLOW_OVERPRESSURE_M3H
        pp = constants.AIRFLOW_PER_PERSON_M3H
        pk = constants.PEOPLE_PEAK
        _patch("AIRFLOW_MIN_M3H", ov)
        _patch("AIRFLOW_MAX_M3H", (ov + pk * pp) * 1.10)
        _patch("P_FAN_RATED_W",
               constants.AIRFLOW_MAX_M3H / 3600.0 * constants.DP_AHU_PA / constants.ETA_FAN)
        _patch("T_HEAT_HIGH_C", constants.T_HEAT_LOW_C + 2.0)

    # --- Sidebar sliders ---
    st.sidebar.header("Parameters")
    h    = st.sidebar.select_slider("Heating setpoint (°C)",     options=[14.0, 16.0, 18.0],       value=16.0)
    a    = st.sidebar.select_slider("Airflow per person (m³/h)", options=[18.0, 25.0, 30.0],       value=18.0)
    hw   = st.sidebar.select_slider("Hot water max (°C)",        options=[40.0, 45.0, 50.0, 55.0], value=40.0)
    s    = st.sidebar.select_slider("Stair cold threshold (°C)", options=[3.0, 5.0, 7.0, 10.0],    value=7.0)
    blow = st.sidebar.select_slider("Blow heat (°C)",            options=[28.0, 30.0, 32.0],       value=28.0)

    # --- Run on button press ---
    if st.button("Run simulation"):
        _patch("T_HEAT_LOW_C",           h)
        _patch("AIRFLOW_PER_PERSON_M3H", a)
        _patch("T_HW_SUPPLY_MAX",        hw)
        _patch("T_STAIR_COLD_C",         s)
        _patch("T_BLOW_HEAT_C",          blow)
        _patch_derived()

        with st.spinner("Running full-year simulation..."):
            r = run_simulation()

        em = r["em"]
        c  = r["comfort"]

        # --- Metric cards ---
        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Annual Cost", f"{em['cost_annual_eur']:.0f} €")
        col2.metric("Energy",      f"{em['E_annual_kWh']:.0f} kWh")
        col3.metric("CO₂",         f"{em['CO2_annual_kgCO2']:.0f} kg")
        col4.metric("Discomfort",  f"{c['combined_discomfort_pct']:.1f} %")

        # --- Energy breakdown bar ---
        fig = go.Figure([go.Bar(
            x=["Heating", "Cooling", "Fans", "Curtain"],
            y=[em["E_heat_total_kWh"], em["E_cool_total_kWh"],
               em["E_fan_total_kWh"],  em["E_curtain_total_kWh"]],
            marker_color=["tomato", "steelblue", "mediumseagreen", "orange"],
        )])
        fig.update_layout(yaxis_title="Energy (kWh)", showlegend=False)
        st.plotly_chart(fig, width="stretch")


# =============================================================================
# PAGE 2 — SOBOL SENSITIVITY
# Variance-based ranking of 15 parameters across 3 metrics.
# =============================================================================

elif page == "Sobol Sensitivity":

    st.subheader("Sobol Sensitivity Analysis")

    @st.cache_data
    def load_sobol():
        import numpy as np
        from SALib.analyze import sobol

        data = np.load("data/processed/sobol_results.npz", allow_pickle=True)

        PROBLEM = {
            "num_vars": 15,
            "names": [
                "T_HEAT_LOW_C", "T_COOL_FIXED_C", "T_HW_SUPPLY_MAX",
                "AIRFLOW_PER_PERSON_M3H", "T_DEAD_LOW_C", "T_DEAD_HIGH_C",
                "T_BLOW_HEAT_C", "T_BLOW_COOL_C", "FRAC_RETURN_AIR",
                "T_NIGHT_SETBACK_C", "T_STAIR_COLD_C", "T_HW_EXT_HIGH_C",
                "T_HW_EXT_HYST_C", "T_CW_EXT_LOW_C", "T_CW_EXT_HYST_C",
            ],
            "bounds": [
                [16, 22], [24, 28], [40, 55], [15, 35], [12, 18], [20, 25],
                [28, 35], [13, 17], [0.50, 0.85], [3, 8], [3, 10], [13, 18],
                [10, 15], [24, 28], [25, 29],
            ],
        }

        Si_cost    = sobol.analyze(PROBLEM, data["Y_cost"],    calc_second_order=True)
        Si_comfort = sobol.analyze(PROBLEM, data["Y_comfort"], calc_second_order=True)
        Si_co2     = sobol.analyze(PROBLEM, data["Y_co2"],     calc_second_order=True)

        return PROBLEM["names"], Si_cost, Si_comfort, Si_co2

    names, Si_cost, Si_comfort, Si_co2 = load_sobol()

    metric = st.selectbox(
        "Metric",
        ["Annual Cost (€)", "Combined Discomfort (%)", "CO₂ (kg)"],
    )

    Si = {
        "Annual Cost (€)":         Si_cost,
        "Combined Discomfort (%)": Si_comfort,
        "CO₂ (kg)":                Si_co2,
    }[metric]

    fig = go.Figure([
        go.Bar(name="S1 (direct effect)", x=names, y=Si["S1"], marker_color="steelblue"),
        go.Bar(name="ST (total effect)",  x=names, y=Si["ST"], marker_color="tomato"),
    ])
    fig.update_layout(
        barmode="group",
        yaxis_title="Sobol Index",
        xaxis_tickangle=-45,
        yaxis=dict(range=[0, 1]),
    )
    st.plotly_chart(fig, width="stretch")


# =============================================================================
# PAGE 3 — PARETO EXPLORER
# Browse precomputed grid: compare a selected config against the optimal one.
# =============================================================================

elif page == "Pareto Explorer":

    st.subheader("Pareto Explorer")

    # --- Load precomputed grid (cached) ---
    @st.cache_data
    def load_pareto():
        return pd.read_csv("data/processed/pareto_all_configs.csv")

    df = load_pareto()

    # --- Optimal config (fixed reference) ---
    OPTIMAL = {
        "T_HEAT_LOW_C":           16,
        "AIRFLOW_PER_PERSON_M3H": 18,
        "T_HW_SUPPLY_MAX":        40,
        "T_STAIR_COLD_C":         7,
        "T_BLOW_HEAT_C":          28,
    }

    optimal_row = df[
        (df["T_HEAT_LOW_C"]           == OPTIMAL["T_HEAT_LOW_C"]) &
        (df["AIRFLOW_PER_PERSON_M3H"] == OPTIMAL["AIRFLOW_PER_PERSON_M3H"]) &
        (df["T_HW_SUPPLY_MAX"]        == OPTIMAL["T_HW_SUPPLY_MAX"]) &
        (df["T_STAIR_COLD_C"]         == OPTIMAL["T_STAIR_COLD_C"]) &
        (df["T_BLOW_HEAT_C"]          == OPTIMAL["T_BLOW_HEAT_C"])
    ]

    # --- Sidebar sliders (snap to grid values via select_slider) ---
    st.sidebar.header("Configuration")
    h = st.sidebar.select_slider(
        "Heating setpoint (°C)",
        options=sorted(df["T_HEAT_LOW_C"].unique()),
        value=OPTIMAL["T_HEAT_LOW_C"],
    )
    a = st.sidebar.select_slider(
        "Airflow per person (m³/h)",
        options=sorted(df["AIRFLOW_PER_PERSON_M3H"].unique()),
        value=OPTIMAL["AIRFLOW_PER_PERSON_M3H"],
    )
    hw = st.sidebar.select_slider(
        "Hot water supply max (°C)",
        options=sorted(df["T_HW_SUPPLY_MAX"].unique()),
        value=OPTIMAL["T_HW_SUPPLY_MAX"],
    )
    s = st.sidebar.select_slider(
        "Staircase cold threshold (°C)",
        options=sorted(df["T_STAIR_COLD_C"].unique()),
        value=OPTIMAL["T_STAIR_COLD_C"],
    )
    blow = st.sidebar.select_slider(
        "Blow heat setpoint (°C)",
        options=sorted(df["T_BLOW_HEAT_C"].unique()),
        value=OPTIMAL["T_BLOW_HEAT_C"],
    )

    # --- Look up the row matching the slider values ---
    filtered = df[
        (df["T_HEAT_LOW_C"]           == h) &
        (df["AIRFLOW_PER_PERSON_M3H"] == a) &
        (df["T_HW_SUPPLY_MAX"]        == hw) &
        (df["T_STAIR_COLD_C"]         == s) &
        (df["T_BLOW_HEAT_C"]          == blow)
    ]

    # --- Display ---
    if filtered.empty:
        st.warning("No matching configuration found.")
    else:
        sel = filtered.iloc[0]      # selected config
        opt = optimal_row.iloc[0]   # optimal config

        # --- Top row: Pareto scatter + CO₂ comparison ---
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Pareto Front")

            # Full cloud of all configs in grey
            fig = px.scatter(
                df,
                x="cost_eur",
                y="combined_discomfort_pct",
                color_discrete_sequence=["lightgrey"],
                labels={
                    "cost_eur": "Annual Cost (€)",
                    "combined_discomfort_pct": "Discomfort (%)",
                },
            )
            fig.update_traces(marker=dict(size=4), showlegend=False)

            # Optimal — green star
            fig.add_trace(go.Scatter(
                x=[opt["cost_eur"]],
                y=[opt["combined_discomfort_pct"]],
                mode="markers",
                marker=dict(color="green", size=12, symbol="star"),
                name="Optimal",
            ))

            # Selected — red dot
            fig.add_trace(go.Scatter(
                x=[sel["cost_eur"]],
                y=[sel["combined_discomfort_pct"]],
                mode="markers",
                marker=dict(color="red", size=12, symbol="circle"),
                name="Selected",
            ))

            st.plotly_chart(fig, width="stretch")

        with col2:
            st.subheader("CO₂ Emissions")
            fig2 = go.Figure([
                go.Bar(name="Optimal",  x=["CO₂ (kg)"], y=[opt["CO2_kg"]], marker_color="green"),
                go.Bar(name="Selected", x=["CO₂ (kg)"], y=[sel["CO2_kg"]], marker_color="red"),
            ])
            fig2.update_layout(barmode="group")
            st.plotly_chart(fig2, width="stretch")

        # --- Bottom: energy breakdown (optimal vs selected) ---
        st.subheader("Energy Breakdown")
        categories = ["Heating", "Cooling", "Fans", "Curtain"]
        opt_vals = [opt["E_heat_kWh"], opt["E_cool_kWh"], opt["E_fan_kWh"], opt["E_curtain_kWh"]]
        sel_vals = [sel["E_heat_kWh"], sel["E_cool_kWh"], sel["E_fan_kWh"], sel["E_curtain_kWh"]]

        fig3 = go.Figure([
            go.Bar(name="Optimal",  x=categories, y=opt_vals, marker_color="green"),
            go.Bar(name="Selected", x=categories, y=sel_vals, marker_color="red"),
        ])
        fig3.update_layout(barmode="group", yaxis_title="Energy (kWh)")
        st.plotly_chart(fig3, width="stretch")