import matplotlib.pyplot as plt
from geothermal_orc_design import multivariate_optimization
from collections import OrderedDict
from znes_plotting import plot, shared
import pandas as pd
import numpy as np
import itertools

import json
import sys


cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/test_multi.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = multivariate_optimization(**input_data)

fluid=list(result.keys())

for a in range(len(fluid)):
    df = pd.DataFrame(list(result.values())[a])
    df['net power output'] = -df['net power output']/1e6
    # Plot the surface.
    X=df['T_before_tur']
    Y=df['IHE_sizing']
    Z=df['dT_air']
    C=df['net power output']
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')  
    p = ax.scatter(X, Y, Z, c=C, cmap=plt.cm.jet)
    ax.set_xlabel('Turbine inlet emperature (°C)')
    ax.set_ylabel('IHE sizing factor (-)')
    ax.set_zlabel('$\Delta T_{air}$ (°C)')
    ax.yaxis.label.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.zaxis.label.set_size(18)
    ax.tick_params(axis="both", labelsize=15)
    cbar = fig.colorbar(p, ax=ax)
    cbar.set_label("Net power output (MW)", size=18)
    cbar.ax.tick_params(labelsize=15)
#    plt.show()
    plt.savefig('Multivariate_optimization_' + fluid[a] + '.pdf')


#for fluid, data in result.items():
#    combinations = list(itertools.combinations(input_data['variables'].keys(), 2))
#
#    for combination in combinations:
#
#        ax = plot.scatter(
#            data=data[data['valid']], x=combination[0], y=combination[1],
#            xlabel=input_data['variables'][combination[0]]['label'],
#            ylabel=input_data['variables'][combination[1]]['label'],
#            colormap=plt.cm.get_cmap('RdYlBu'), c=input_data['objective'])
#
#    # Create multipage PDF from plots
#    shared.create_multipage_pdf(fluid + '_multivariate_scatterplots.pdf')
