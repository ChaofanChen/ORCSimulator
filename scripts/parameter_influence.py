from geothermal_orc import single_parameter_influence
import json
import sys

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

#%% Turbine inlet temperature without IHE

with open(cur_dir + '/T_tur_influence.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)

#%% IHE sizing factor at optimal turbine inlet temperature

for fluid in ['R245fa', 'R600', 'R245ca', 'Isopentane']:

    with open(cur_dir + '/IHE_sizing_influence_' + fluid + '.json', 'r') as f:
        input_data = json.load(f)
        f.close()

    result = single_parameter_influence(**input_data)

#%% Exemplary Ts plots at maximum IHE size and optimal turbine inlet temperature

for fluid in ['R245ca', 'R600']:

    with open(cur_dir + '/' + fluid + '_Ts_plot.json', 'r') as f:
        input_data = json.load(f)
        f.close()

    result = single_parameter_influence(**input_data)

#%% Condenser air temperature change

with open(cur_dir + '/T_cond_influence.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)

#%% Variation of geosteam share

with open(cur_dir + '/Low_geo_steam.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)

#%% Variation of turbine inlet temperature with IHE size to realise 70 °C re-injection

with open(cur_dir + '/IHE_install.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)
