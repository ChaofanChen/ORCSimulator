import matplotlib.pyplot as plt
from geothermal_orc_design import single_optimization, single_parameter_influence, multivariate_optimization
from collections import OrderedDict

import itertools

import json
import sys


cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/test.json', 'r') as f:
    input_data = json.load(f)
    f.close()

single_parameter_result = single_parameter_influence(**input_data)
single_parameter_opt_result = single_optimization(**input_data)

data = input_data.copy()
del data['boundary_conditions']['IHE_sizing']
single_parameter_influence_at_opt = {}
for fluid in input_data['working_fluid_list']:
    data['working_fluid_list'] = [fluid]
    data['boundary_conditions']['T_before_tur'] = single_parameter_opt_result['T_before_tur'].loc[fluid, 'T_1']
    data['variables'] = {
        'IHE_sizing': {
            "max": 1,
            "min": 0,
            "tol": 1e-2,
            "unit": "-",
            "label": "Internal heat exchanger sizing factor"
        }
    }

    single_parameter_influence_at_opt[fluid] = single_parameter_influence(**data)[fluid]

data['working_fluid_list'] = input_data['working_fluid_list']
data['variables']['T_before_tur'] = input_data['variables']['T_before_tur']
data['variables']['dT_air'] = {
    "max": 35,
    "min": 10,
    "tol": 1e-2,
    "unit": "K",
    "label": "Condenser air temperature increase"
}

del data['boundary_conditions']['T_before_tur']
del data['boundary_conditions']['dT_air']

data['objective'] = 'net power output'

# multivariate_opt_result = multivariate_optimization(**data)

print(single_parameter_result)
print(single_parameter_opt_result)
print(single_parameter_influence_at_opt)
# print(multivariate_opt_result)
