import matplotlib.pyplot as plt
from geothermal_orc_design import multivariate_optimization
from collections import OrderedDict
from znes_plotting import plot, shared
import pandas as pd
import numpy as np
import itertools
from matplotlib.ticker import MaxNLocator
import json
import sys


cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/test_multi.json', 'r') as f:
    input_data = json.load(f)
    f.close()

opt_results = pd.DataFrame(columns=['fluid', 'net_power', 'T_before_tur', 'dT_air'])

result, opt_results = multivariate_optimization(**input_data)

fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
ax.bar(opt_results['fluid'], opt_results['net_power'], color='blue', width=0.4)#, marker="o", linewidth=2
for i, v in enumerate(np.array(opt_results['net_power'])):
    ax.text(i - 0.25, v + 0.1, str(round(v, 2)), color='black', fontweight='bold', fontsize=12)
ax.set(xlabel= 'Working fluid', ylabel='Net power output (MW)')
ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18)
plt.ylim(10, 13)
ax.tick_params(axis="x", labelsize=15, width=0.5)
ax.tick_params(axis="y", labelsize=15, width=0.5)
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(1)
fig.autofmt_xdate()
ax.grid(b=True, which='major', linewidth=0.5)
plt.savefig('diff_working_fluid_net_power.pdf')
plt.show()

print(opt_results.to_latex(escape=False, na_rep='-', float_format='%.2f'))


#fluid=list(result.keys())
#
#for a in range(len(fluid)):
#    df = pd.DataFrame(list(result.values())[a])
#
#    # Plot the surface.
#    X=list(df['T_before_tur'])
#    Y=list(df['IHE_sizing'])
#    Z=list(df['dT_air'])
#    C=list(df['net power output'])
#
##    from mpl_toolkits.mplot3d import Axes3D
##    fig = plt.figure()
##    ax = Axes3D(fig)
##    surf = ax.plot_trisurf(X, Z, Y, cmap=plt.cm.jet)
##    surf = ax.contour3D(X, Z, C, cmap=plt.cm.jet, linewidth=0.01)
##    fig.colorbar(surf, shrink=0.5, aspect=5)
##    plt.show()
#
#    fig = plt.figure(figsize=(16, 12), dpi=100)
#    ax = fig.add_subplot(111, projection='3d')
#    surf = ax.scatter(X, Z, C, c=C, cmap=plt.cm.jet, s=80)
#    ax.set_xlabel('Turbine inlet temperature (°C)')
#    ax.set_ylabel('$\Delta T_{air}$ (°C)')
#    ax.set_zlabel('Net power output (MW)')
#    ax.xaxis.labelpad=12
#    ax.yaxis.labelpad=12
#    ax.zaxis.labelpad=12
#
#    ax.yaxis.label.set_size(20)
#    ax.xaxis.label.set_size(20)
#    ax.zaxis.label.set_size(20)
#    ax.tick_params(axis="both", labelsize=18)
#    cbar = fig.colorbar(surf, ax=ax)
#    cbar.set_label("Net power output (MW)", size=20)
#    cbar.ax.tick_params(labelsize=18)
#    plt.savefig('Multivariate_optimization_' + fluid[a] + '.pdf')
#    plt.show()


#%%
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
