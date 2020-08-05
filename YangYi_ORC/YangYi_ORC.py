from tespy.connections import connection, bus
from tespy.tools import char_line
from tespy.networks import network
from tespy.components import heat_exchanger, pump, turbine, source, sink, cycle_closer, splitter, merge, condenser
from tespy.components.customs import orc_evaporator
from CoolProp.CoolProp import PropsSI
from CoolProp.Plots import SimpleCompressionCycle
import matplotlib.pyplot as plt
from CoolProp.Plots import PropertyPlot
import CoolProp
import CoolProp as CP
import pandas as pd
import numpy as np
from tespy.tools import logger
import logging
mypath = logger.define_logging(
log_path=True, log_version=True, timed_rotating={'backupCount': 4},
screen_level=logging.WARNING, screen_datefmt = "no_date")
print('CoolProp ver:%s'%(CoolProp.__version__))
# define basic cycle
fluids = ['water', 'Isopentane', 'Air']
nw = network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
# input parameters (the mass flow rate of cooling air should be adjusted
# based on the temperature of the geo-fluid for stable calculation)
# geo-fluid part
# mass_flow_rate_brine = 190.56 # kg/s
volume_flow_rate_brine = 700.23 # m3/h
mass_flow_rate_steam = 20.28
T_brine_in = 144.8
T_reinjection = 70.8
# cooling air part
# mass_flow_rate_air = 6284.6 # 6241.5
T_air = 0.5
p_air = 0.61
# calculation secondary variables
p_before_turbine = PropsSI('P', 'T', T_brine_in+273.15-22.6, 'Q', 1, 'Isopentane')/1e5
rho_brine_in=PropsSI('D', 'T', T_brine_in+273.15, 'Q', 0, 'water')
mass_flow_rate_brine = volume_flow_rate_brine*rho_brine_in /3600

p_steam_in = PropsSI('P', 'T', T_brine_in+273.15, 'Q', 1, 'water')/1e5

#T=PropsSI('T', 'P', 0.8e5, 'Q', 0, 'Isopentane')-273.15
# main components
evaporator = orc_evaporator('orc_evaporator')
pump_c = pump('condensate pump')
merge = merge('geo-fluid merge point')
preheater = heat_exchanger('preheater')
turbine = turbine('turbine')
ihe = heat_exchanger('internal heat exchanger')
condenser = condenser('condenser')
pump = pump('feeding pump')
# working fluid
source_wf = source('working fluid source')
sink_wf = sink('working fluid sink')
close_cycle = cycle_closer('cycle closer before preheater')
#brine
source_s = source('steam source')
source_b = source('brine source')
sink_s = sink('steam sink')
sink_b = sink('brine sink')
# cooling air
source_ca = source('cooling air source')
sink_ca = sink('cooling air sink')

# connections
# main cycle
preheater_wf_in = connection(close_cycle, 'out1', preheater, 'in2')
preheater_evaporator = connection(preheater, 'out2', evaporator, 'in3')
evaporator_turbine = connection(evaporator, 'out3', turbine, 'in1')
turbine_ihe = connection(turbine, 'out1', ihe, 'in1')
ihe_condenser = connection(ihe, 'out1', condenser, 'in1')
condenser_pump = connection(condenser, 'out1', pump, 'in1')
pump_ihe = connection(pump, 'out1', ihe, 'in2')
ihe_wf_out = connection(ihe, 'out2', close_cycle, 'in1')
nw.add_conns(preheater_wf_in, preheater_evaporator, evaporator_turbine, turbine_ihe, ihe_condenser, condenser_pump, pump_ihe, ihe_wf_out)
# geo-steam cycle
evaporator_steam_in = connection(source_s, 'out1', evaporator, 'in1')
evaporator_pump = connection(evaporator, 'out1', pump_c, 'in1')
pump_sink_s = connection(pump_c, 'out1', merge, 'in1')
# geo-brine cycle
evaporator_brine_in = connection(source_b, 'out1', evaporator, 'in2')
evaporator_sink_b = connection(evaporator, 'out2', merge, 'in2')
merge_preheater = connection(merge, 'out1', preheater, 'in1')
preheater_sink = connection(preheater, 'out1', sink_b, 'in1')
nw.add_conns(evaporator_steam_in, evaporator_pump, pump_sink_s, evaporator_brine_in, evaporator_sink_b, merge_preheater, preheater_sink)
# cooling air cycle
ca_in = connection(source_ca, 'out1', condenser, 'in2')
ca_out = connection(condenser, 'out2', sink_ca, 'in1')
nw.add_conns(ca_in, ca_out)

# parametrization of components
evaporator.set_attr(pr1=0.81181818, pr2=0.970588, pr3=1)
preheater.set_attr(pr1=0.949494, pr2=0.955752)
pump_c.set_attr(pr=2.4480712, eta_s=0.8)
turbine.set_attr(pr=0.114012184, eta_s=0.85, design=['eta_s', 'pr'])
pump.set_attr(eta_s=0.9)
ihe.set_attr(pr1=0.849056603, pr2=0.957627118)
condenser.set_attr(pr1=0.8889, pr2=1)

# busses
# characteristic function for generator efficiency
# x = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
# y = np.array([0, 0.88, 0.89, 0.90, 0.91, 0.976, 0.91])
# gen = char_line(x=x, y=y)
# # motor of pump has a constant efficiency
# power = bus('total output power')
# power.add_comps({'c': turbine, 'p': 'P', 'char': gen})
# nw.add_busses(power)

# parametrization of connections
preheater_evaporator.set_attr(p=p_before_turbine, fluid={'water': 0, 'Isopentane': 1, 'Air': 0})

evaporator_steam_in.set_attr(T=T_brine_in, p=p_steam_in, m=mass_flow_rate_steam, state='g', fluid={'water': 1, 'Isopentane': 0, 'Air':0})
evaporator_brine_in.set_attr(T=T_brine_in, m=mass_flow_rate_brine, fluid={'water': 1, 'Isopentane': 0, 'Air':0})
preheater_sink.set_attr(T=T_reinjection)
# evaporator_pump.set_attr(Td_bp=-5)
# evaporator_sink_b.set_attr(T=T_brine_in-22)

# air cooling connections
ca_in.set_attr(T=T_air, p=p_air, fluid={'water': 0, 'Isopentane': 0, 'Air': 1})
# ca_out.set_attr(T=T_air + 15)

# parametrization of components
# The parameter "ttd_u" will
# constrain the temperature difference between
# the hot inlet and the cold outlet.
# This value is based on experience from other
# calculations. It will only influence the heat
# ejected by the air condenser.
ihe.set_attr(ttd_u=23.2)
# Here the hot inlet and cold outlet temperature
# difference is constrained. The value is based on
# Yangyi monitoring data.
preheater.set_attr(ttd_u=5.1)
condenser.set_attr(ttd_u=6.98)
# solving
mode = 'design'
save_path = 'yangyi'
nw.solve(mode=mode, init_path=save_path)
nw.print_results()
##%-----------------------------------------------------------------------------------------------------------
s_before_turbine = PropsSI('S', 'T', evaporator_turbine.T.val + 273.15, 'Q', 1, 'Isopentane')
T_before_turbine = evaporator_turbine.T.val

s_after_turbine = PropsSI('S', 'T', turbine_ihe.T.val + 273.15, 'P', turbine_ihe.p.val * 100000, 'Isopentane')
T_after_turbine = turbine_ihe.T.val

s_before_condenser = PropsSI('S', 'T', ihe_condenser.T.val + 273.15, 'P', ihe_condenser.p.val * 100000, 'Isopentane')
T_before_condenser = ihe_condenser.T.val

s_after_condenser = PropsSI('S', 'T', condenser_pump.T.val + 273.15, 'Q', 0, 'Isopentane')
T_after_condenser = condenser_pump.T.val

s_after_pump = PropsSI('S', 'T', pump_ihe.T.val + 273.15, 'P', pump_ihe.p.val * 100000, 'Isopentane')
T_after_pump = pump_ihe.T.val

s_after_ihe = PropsSI('S', 'T', ihe_wf_out.T.val + 273.15, 'P', ihe_wf_out.p.val * 100000, 'Isopentane')
T_after_ihe = ihe_wf_out.T.val

s_after_preheater = PropsSI('S', 'T', preheater_evaporator.T.val + 273.15, 'P', preheater_evaporator.p.val * 100000, 'Isopentane')
T_after_preheater = preheater_evaporator.T.val

print('Power_output (MW):', turbine.P.val / 1e6)

state = CP.AbstractState('HEOS', 'Isopentane')
T_crit = state.trivial_keyed_output(CP.iT_critical)
df = pd.DataFrame(columns=['s_l', 's_g', 's_iso_P0', 's_iso_P1', 's_iso_P_top', 's_iso_P_bottom'])
P0 = condenser_pump.p.val * 100000
P1 = p_before_turbine * 100000
T_range = np.geomspace(273.15, T_crit, 1000)
for T in T_range:
    df.loc[T, 's_l'] = PropsSI('S', 'T', T, 'Q', 0, 'Isopentane')
    df.loc[T, 's_g'] = PropsSI('S', 'T', T, 'Q', 1, 'Isopentane')
    df.loc[T, 's_iso_P0'] = PropsSI('S', 'T', T, 'P', P0, 'Isopentane')
    df.loc[T, 's_iso_P1'] = PropsSI('S', 'T', T, 'P', P1, 'Isopentane')

T_range_evaporator = np.geomspace(preheater_evaporator.T.val + 273.15, evaporator_turbine.T.val + 273.15+0.1, 100)
for T in T_range_evaporator:
    df.loc[T, 's_iso_P_top'] = PropsSI('S', 'T', T, 'P', P1, 'Isopentane')

T_steam_wf_low_P = PropsSI('T', 'P', P0, 'Q', 1, 'Isopentane')
s_steam_wf_low_P = PropsSI('S', 'P', P0, 'Q', 1, 'Isopentane')
T_range_condenser = np.geomspace(T_steam_wf_low_P + 0.1, condenser_pump.T.val + 273.15-0.1, 100)
for T in T_range_condenser:
    df.loc[T, 's_iso_P_bottom'] = PropsSI('S', 'T', T, 'P', P0, 'Isopentane')
# print(df)

fig, ax = plt.subplots()
ax.plot(df['s_g'], df.index - 273.15, color='black')
ax.plot(df['s_l'], df.index - 273.15, color='black')
ax.plot(df['s_iso_P0'], df.index - 273.15, color='green')
ax.plot(df['s_iso_P1'], df.index - 273.15, color='green')
ax.plot(df['s_iso_P_top'], df.index - 273.15, color='red')
ax.plot(df['s_iso_P_bottom'], df.index - 273.15, color='red')

Temp = [T_before_turbine, T_after_turbine, T_before_condenser, T_steam_wf_low_P-273.15, T_after_condenser, T_after_ihe, T_after_preheater]
entropy = [s_before_turbine, s_after_turbine, s_before_condenser, s_steam_wf_low_P, s_after_condenser, s_after_ihe, s_after_preheater]
n = ['before turbine', 'after turbine', 'before condenser', ' ', 'after condenser', 'after ihe', 'after preheater']

ax.scatter(entropy, Temp, color='red')
for i, txt in enumerate(n):
    ax.annotate(txt, (entropy[i], Temp[i])) #, ha='center'
for i in range (0, 3, 1):
    plt.plot(entropy[i:i+2], Temp[i:i+2], 'ro-', lw=2)
for i in range (4, 6, 1):
    plt.plot(entropy[i:i+2], Temp[i:i+2], 'ro-', lw=2)

ax.set(xlabel='Specific entropy [J/kg K]', ylabel='Temperature [dC]',
       title='T,s Graph for working fluid')
ax.grid()
plt.savefig('ts_plot_with_ORC_cycle.png')
plt.show()