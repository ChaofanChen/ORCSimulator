import matplotlib.pyplot as plt
from geothermal_orc_design import single_parameter_influence
from collections import OrderedDict
import pandas as pd
import numpy as np
import CoolProp as CP
from CoolProp.CoolProp import PropsSI as PSI

import json
import sys

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/diff_fluid_performance.json', 'r') as f:
    input_data = json.load(f)
    f.close()


result = single_parameter_influence(**input_data)

### plot figures
fluid=list(result.keys())

fig, ax = plt.subplots(figsize=(8, 6), dpi=100)

for a in range(len(fluid)):
    df = pd.DataFrame(list(result.values())[a])
    df['gross power output'] = -df['gross power output']/1e6
    df['net power output'] = -df['net power output']/1e6


    ax.plot(fluid[a], df['gross power output'], marker='x', label=fluid[a])

    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.tick_params(axis='y', colors='blue')

    for i, v in enumerate(np.array(df['gross power output'])):
        ax.text(i+a+0.01, v+0.01, str(round(v, 2)), color='black',
                fontweight='bold', fontsize=12)
ax.set(xlabel= 'Fluid', ylabel='Gross power output (MW)')
ax.grid()
plt.ylim(2.8, 3.1)
ax.legend()
plt.savefig('diff_fluid_performance.png')
