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
#    plt.ylim(40, 100)
#    plt.xlim(50, 130)
    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_size(18)
    ax2.tick_params(axis="y", labelsize=15)
    ax.grid()
#    plt.savefig('T_tur_with_' + fluid[a] + '.pdf')

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
    ax.plot(df['T_1'], df['gross power output'], color='blue', marker="o", label='Gross power output')
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
    ax.legend(frameon=False, prop={'size': 12})
    plt.savefig('Low_geo_steam_' + fluid[a] + '.pdf')
    df.to_csv('Low_geo_steam_' + fluid[a] + '.csv')

#%%

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/R245ca_Ts_plot.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)

print(result)

fluid=list(result.keys())[0]
df = pd.DataFrame(list(result.values())[0])

T_before_turbine = df.iloc[0]['T_1']
s_before_turbine = df.iloc[0]['s_1']

T_after_turbine = df.iloc[0]['T_2']
s_after_turbine = df.iloc[0]['s_2']

T_before_condenser = df.iloc[0]['T_3']
s_before_condenser = df.iloc[0]['s_3']
T_after_condenser = df.iloc[0]['T_4']
s_after_condenser = df.iloc[0]['s_4']

T_after_pump = df.iloc[0]['T_5']
s_after_pump = df.iloc[0]['s_5']

T_after_ihe = df.iloc[0]['T_6']
s_after_ihe = df.iloc[0]['s_6']

T_after_preheater = df.iloc[0]['T_7']
s_after_preheater = df.iloc[0]['s_7']

state = CP.AbstractState('HEOS', fluid)
T_crit = state.trivial_keyed_output(CP.iT_critical)
df_data = pd.DataFrame(columns=['s_l', 's_g', 's_iso_P0', 's_iso_P1', 's_iso_P_top', 's_iso_P_bottom'])

P0 = PSI('P', 'T', T_after_condenser + 273.15, 'Q', 0, fluid)
P1 = PSI('P', 'T', T_before_turbine + 273.15, 'Q', 1, fluid)

T_range = np.linspace(273.15, T_crit, 1000)
for T in T_range:
    df_data.loc[T, 's_l'] = PSI('S', 'T', T, 'Q', 0, fluid)
    df_data.loc[T, 's_g'] = PSI('S', 'T', T, 'Q', 1, fluid)
    df_data.loc[T, 's_iso_P0'] = PSI('S', 'T', T, 'P', P0, fluid)
    df_data.loc[T, 's_iso_P1'] = PSI('S', 'T', T, 'P', P1, fluid)

T_range_evaporator = np.linspace(T_after_preheater + 273.15, T_before_turbine + 273.15 + 0.5, 10)
for T in T_range_evaporator:
    df_data.loc[T, 's_iso_P_top'] = PSI('S', 'T', T, 'P', P1, fluid)

T_range_condenser = np.linspace(T_before_condenser + 273.15, T_after_condenser + 273.15 - 0.5, 20)
for T in T_range_condenser:
    df_data.loc[T, 's_iso_P_bottom'] = PSI('S', 'T', T, 'P', P0, fluid)
#print(df)

fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
ax.plot(df_data['s_g'], df_data.index - 273.15, color='black')
ax.plot(df_data['s_l'], df_data.index - 273.15, color='black')
ax.plot(df_data['s_iso_P0'], df_data.index - 273.15, color='green')
ax.plot(df_data['s_iso_P1'], df_data.index - 273.15, color='green')
ax.plot(df_data['s_iso_P_top'], df_data.index - 273.15, color='red')
ax.plot(df_data['s_iso_P_bottom'], df_data.index - 273.15, color='red')

Temp = [T_before_turbine, T_after_turbine, T_before_condenser, T_after_condenser, T_after_pump,
        T_after_ihe, T_after_preheater]
entropy = [s_before_turbine, s_after_turbine, s_before_condenser, s_after_condenser, s_after_pump, s_after_ihe,
           s_after_preheater]
n = ['1', '2', '3', ' ', '5', '6', '7']

ax.scatter(entropy, Temp, color='red', s=0.5)
for i, txt in enumerate(n):
    ax.annotate(txt, (entropy[i], Temp[i]), fontsize=16, textcoords="offset points", xytext=(0,6), horizontalalignment='right')  # , ha='center'
for i in range(0, 2, 1):
    plt.plot(entropy[i:i + 2], Temp[i:i + 2], 'ro-', lw=2)
for i in range(4, 6, 1):
    plt.plot(entropy[i:i + 2], Temp[i:i + 2], 'ro-', lw=2)

# geo-source
ax.scatter(df['s_30'] * 0.1 + df['s_32'] * 0.9 - 250, df['T_30']+10)
ax.scatter(df['s_35'], df['T_35'] + 20)
# cooling air
ax.scatter(df['s_20']-3000, df['T_20'])
ax.scatter(df['s_22']-2100, df['T_22'])

plt.ylim(0, 180)
ax.set(xlabel='Specific entropy', ylabel='Temperature')
#ax.set(xlabel='Specific entropy (J/kg K)', ylabel='Temperature (°C)')
ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18)
#ax.tick_params(axis="x", labelsize=15)
#ax.tick_params(axis="y", labelsize=15)
#ax.grid()
# Turn off tick labels
ax.set_yticklabels([])
ax.set_xticklabels([])
#plt.savefig('ORC_Ts_plot_' + fluid + '.pdf')
plt.savefig('ORC_Ts_plot.pdf')
plt.show()

#%%

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/IHE_install.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_parameter_influence(**input_data)

### plot figures
fluid=list(result.keys())

for a in range(len(fluid)):
    df = pd.DataFrame(list(result.values())[a])
    df['gross power output'] = -df['gross power output']/1e6
    df['net power output'] = -df['net power output']/1e6
    df['Q_internal heat exchanger'] = -df['Q_internal heat exchanger']/1e6

    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    ax.plot(df['T_1'], df['net power output'], color='blue', marker="v", label='Net power output')
    ax.set(xlabel= 'Turbine inlet temperature with ' + fluid[a] + ' (°C)', ylabel='Net power output (MW)')
    plt.ylim(12, 13)
    ax2=ax.twinx()
    ax2.plot(df['T_1'], df['Q_internal heat exchanger'], color='black', marker="*", label='Q_internal heat exchanger')
    ax2.set(ylabel='Required IHE heat exchange rate (MW)')
    plt.ylim(5.5, 8)
    plt.xlim(125, 131)
    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    ax.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_size(18)
    ax2.tick_params(axis="y", labelsize=15)
    ax.grid()
#    ax.legend()
    plt.savefig('IHE_install' + fluid[a] + '.pdf')
    df.to_csv('IHE_install_' + fluid[a] + '.csv')


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
