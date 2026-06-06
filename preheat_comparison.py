from simulation import run_simulation
import regulation 

preheat_state = ["none", "setback_override", "hc"]

results = {}
for pre in preheat_state:
    regulation.PREHEAT_STRATEGY = pre
    results[pre] = run_simulation()
    
r_none = results["none"]
r_setback = results["setback_override"]
r_hc = results["hc"]

def pct_diff(value, baseline):
    return (value - baseline) / baseline * 100 if baseline else 0.0

print ("\nPreheat strategy: none")
print ("Total energy: %0.2f kWh" % r_none["em"]["E_heat_total_kWh"])
print ("Total cost: %0.2f $ " % r_none["em"]["cost_annual_eur"])
print ("Total CO2: %0.2f kg " % r_none["em"]["CO2_annual_kgCO2"])
print("Total discomfort: %0.2f %%" % r_none["comfort"]["combined_discomfort_pct"])

print ("\nPreheat strategy: setback_override")
print ("Total energy: %0.2f kWh (%+0.1f%% vs none) " % (r_setback["em"]["E_heat_total_kWh"], 
    pct_diff(r_setback["em"]["E_heat_total_kWh"], r_none["em"]["E_heat_total_kWh"])))
print ("Total cost: %0.2f $ (%+0.1f%% vs none)" % (r_setback["em"]["cost_annual_eur"], 
    pct_diff(r_setback["em"]["cost_annual_eur"], r_none["em"]["cost_annual_eur"])))
print ("Total CO2: %0.2f kg (%+0.1f%% vs none)" % (r_setback["em"]["CO2_annual_kgCO2"], 
    pct_diff(r_setback["em"]["CO2_annual_kgCO2"], r_none["em"]["CO2_annual_kgCO2"])))
print ("Total discomfort: %0.2f %% (%+0.1f%% vs none)" % (r_setback["comfort"]["combined_discomfort_pct"], 
    pct_diff(r_setback["comfort"]["combined_discomfort_pct"], r_none["comfort"]["combined_discomfort_pct"])))

print ("\nPreheat strategy: hc")    
print ("Total energy: %0.2f kWh (%+0.1f%% vs none)" % (r_hc["em"]["E_heat_total_kWh"], 
    pct_diff(r_hc["em"]["E_heat_total_kWh"], r_none["em"]["E_heat_total_kWh"])))
print ("Total cost: %0.2f $ (%+0.1f%% vs none)" % (r_hc["em"]["cost_annual_eur"], 
    pct_diff(r_hc["em"]["cost_annual_eur"], r_none["em"]["cost_annual_eur"])))
print ("Total CO2: %0.2f kg (%+0.1f%% vs none)" % (r_hc["em"]["CO2_annual_kgCO2"], 
    pct_diff(r_hc["em"]["CO2_annual_kgCO2"], r_none["em"]["CO2_annual_kgCO2"])))
print ("Total discomfort: %0.2f %% (%+0.1f%% vs none)" % (r_hc["comfort"]["combined_discomfort_pct"], 
    pct_diff(r_hc["comfort"]["combined_discomfort_pct"], r_none["comfort"]["combined_discomfort_pct"])))

