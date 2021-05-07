import matplotlib.pyplot as plt
from geothermal_orc_design import single_parameter_influence
from collections import OrderedDict
import pandas as pd
import numpy as np

import json
import sys

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/T_tur_influence.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)

print(result)

### plot figures
fluid=list(result.keys())

for a in range(len(fluid)):
    df = pd.DataFrame(list(result.values())[a])
    df['gross power output'] = -df['gross power output']/1e6
    df['net power output'] = -df['net power output']/1e6

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    ax.plot(df['T_1'], df['gross power output'], color='blue', marker="o")
    ax.set(xlabel= 'Turbine inlet temperature with ' + fluid[a] + ' (°C)', ylabel='Gross power output (MW)')
    plt.ylim(5, 18)
    ax2=ax.twinx()
    ax2.plot(df['T_1'], df['T_35'], color='black', marker="*", label='Re-injection temperature')
    ax2.set(ylabel='Re-injection temperature (°C)')
    plt.ylim(40, 100)
    plt.xlim(50, 130)
    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_size(18)
    ax2.tick_params(axis="y", labelsize=15)
    ax.grid()   
    plt.savefig('T_tur_with_' + fluid[a] + '.pdf')

#%%   IHE sizing factor

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/IHE_sizing_influence.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)

print(result)

### plot figures
fluid=list(result.keys())

for a in range(len(fluid)):
    df = pd.DataFrame(list(result.values())[a])
    df['gross power output'] = -df['gross power output']/1e6
    df['net power output'] = -df['net power output']/1e6
    df['Q_internal heat exchanger'] = -df['Q_internal heat exchanger']/1e6

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    ax.plot(df['Q_internal heat exchanger'], df['net power output'], color='blue', marker="o")
    ax.set(xlabel= 'IHE heat exchange rate with ' + fluid[a] + '(MW)', ylabel='Net power output (MW)')
    plt.ylim(8, 14)
    ax2=ax.twinx()
    ax2.plot(df['Q_internal heat exchanger'], df['T_35'], color='black', marker="*", label='Re-injection temperature')
    ax2.set(ylabel='Re-injection temperature (°C)')
    plt.ylim(40, 100)
    plt.xlim(0, 8)
    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_size(18)
    ax2.tick_params(axis="y", labelsize=15)
    ax2.plot(np.linspace(0, 8, 5), (70,70,70,70,70), dashes=[2, 2], linewidth=3, color='red')
    ax.grid()   
    plt.savefig('IHE_sizing_with_' + fluid[a] + '.pdf')
    df.to_csv('IHE_sizing_with_' + fluid[a] + '.csv')


#%%

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/T_cond_influence.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)

print(result)

## plot figures
fluid=list(result.keys())

for a in range(len(fluid)):
    df = pd.DataFrame(list(result.values())[a])
    df['gross power output'] = -df['gross power output']/1e6
    df['net power output'] = -df['net power output']/1e6
    df['net efficiency'] = df['net efficiency'] * 100

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    ax.plot(df['T_4'], df['net power output'], color='blue', marker="o")
    ax.set(xlabel= 'Condenser outlet temperature with ' + fluid[a] + ' (°C)', ylabel='Net power output (MW)')
    plt.ylim(5, 14)
    ax2=ax.twinx()
    ax2.plot(df['T_4'], df['net efficiency'], color='black', marker="*", label='Net efficiency (%)')
    ax2.set(ylabel='Net efficiency (%)')
    plt.ylim(0, 14)
    plt.xlim(20, 45)
    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_size(18)
    ax2.tick_params(axis="y", labelsize=15)
    ax.grid()
    plt.savefig('T_cond_with_' + fluid[a] + '.pdf')

#%%

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/Low_geo_steam.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)

print(result)

### plot figures
fluid=list(result.keys())

for a in range(len(fluid)):
    df = pd.DataFrame(list(result.values())[a])
    df['gross power output'] = -df['gross power output']/1e6
    df['net power output'] = -df['net power output']/1e6

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    ax.plot(df['T_1'], df['gross power output'], color='r', marker="o", label='Gross power output')
    ax.plot(df['T_1'], df['net power output'], color='blue', marker="v", label='Net power output')
    ax.set(xlabel= 'Turbine inlet temperature with ' + fluid[a] + ' (°C)', ylabel='Power output (MW)')
    plt.ylim(0, 12)
    ax2=ax.twinx()
    ax2.plot(df['T_1'], df['T_35'], color='black', marker="*", label='Re-injection temperature')
    ax2.set(ylabel='Re-injection temperature (°C)')
    plt.ylim(40, 110)
    plt.xlim(50, 130)
    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_size(18)
    ax2.tick_params(axis="y", labelsize=15)
    ax.grid()
    ax.legend()
    plt.savefig('Low_geo_steam_' + fluid[a] + '.pdf')
    df.to_csv('Low_geo_steam_' + fluid[a] + '.csv')




#%%
#cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'
#
#with open(cur_dir + '/Q_BEv_influence.json', 'r') as f:
#    input_data = json.load(f)
#    f.close()
#
#result = single_parameter_influence(**input_data)
#
#print(result)
#df_Q_bev = pd.DataFrame(list(result.values())[0])
#df_Q_bev['gross power output'] = -df_Q_bev['gross power output']/1e6
#df_Q_bev['net power output'] = -df_Q_bev['net power output']/1e6
#
#fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
#ax.plot(df_Q_bev['T_1'], df_Q_bev['gross power output'], color='blue', marker="o")
#ax.set(xlabel= 'Turbine inlet temperature (bar)', ylabel='Gross power output (MW)')
##plt.ylim(10, 18)
#ax2=ax.twinx()
#ax2.plot(df_Q_bev['T_1'], df['T_35'], color='black', marker="*", label='Re-injection temperature')
#ax2.set(ylabel='Re-injection temperature (°C)')
##plt.ylim(60, 75)
##plt.xlim(0, 8)
#ax.yaxis.label.set_color('blue')
#ax.yaxis.label.set_size(15)
#ax.xaxis.label.set_size(15)
#ax.tick_params(axis='y', colors='blue')
#ax2.yaxis.label.set_size(15)
#ax.grid()
