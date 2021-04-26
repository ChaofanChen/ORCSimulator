import matplotlib.pyplot as plt
from geothermal_orc_design import single_optimization
from collections import OrderedDict
import pandas as pd
import itertools

import json
import sys


cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/test.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_optimization(**input_data)

df = pd.DataFrame(list(result.values())[0], index = input_data['working_fluid_list'])

df['gross power output'] = -df['gross power output']/1e6
df['net power output'] = -df['net power output']/1e6
df['thermal efficiency'] = df['thermal efficiency']*100
df['net efficiency'] = df['net efficiency']*100
df['Q_geobrine evaporator'] = -df['Q_geobrine evaporator']/1e6
print(df.to_latex(escape=False, na_rep='-', float_format='%.2f'))