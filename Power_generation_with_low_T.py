from tespy.connections import connection
from tespy.tools import char_line
from tespy.networks import network
from tespy.components import heat_exchanger, pump, turbine, source, sink, cycle_closer, splitter, merge, condenser
from tespy.components.customs import orc_evaporator
from CoolProp.CoolProp import PropsSI
import numpy as np
from tespy.tools import logger
import logging
mypath = logger.define_logging(
log_path=True, log_version=True, timed_rotating={'backupCount': 4},
screen_level=logging.WARNING, screen_datefmt = "no_date")
# define basic cycle
fluids = ['water', 'Isopentane', 'Air']
nw = network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
# input parameters (the mass flow rate of cooling air should be adjusted
# based on the temperature of the geo-fluid for stable calculation)
# geo-fluid part
mass_flow_rate_brine = 200
T_brine_in = 120
T_reinjection = 45
# cooling air part
# mass_flow_rate_air = 6284.6 # 6241.5
T_air = 4.7
p_air = 0.91
# calculation secondary variables
p_before_turbine = PropsSI('P', 'T', T_brine_in+273.15-26.8, 'Q', 1, 'Isopentane')/1e5
#T=PropsSI('T', 'P', 0.8e5, 'Q', 0, 'Isopentane')-273.15
# main components
evaporator = heat_exchanger('orc_evaporator')
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
source_b = source('brine source')
sink_b = sink('brine sink')
# cooling air
source_ca = source('cooling air source')
sink_ca = sink('cooling air sink')

# connections
# main cycle
preheater_wf_in = connection(close_cycle, 'out1', preheater, 'in2')
preheater_evaporator = connection(preheater, 'out2', evaporator, 'in2')
evaporator_turbine = connection(evaporator, 'out2', turbine, 'in1')
turbine_ihe = connection(turbine, 'out1', ihe, 'in1')
ihe_condenser = connection(ihe, 'out1', condenser, 'in1')
condenser_pump = connection(condenser, 'out1', pump, 'in1')
pump_ihe = connection(pump, 'out1', ihe, 'in2')
ihe_wf_out = connection(ihe, 'out2', close_cycle, 'in1')
nw.add_conns(preheater_wf_in, preheater_evaporator, evaporator_turbine, turbine_ihe, ihe_condenser, condenser_pump, pump_ihe, ihe_wf_out)
# geo-brine cycle
evaporator_brine_in = connection(source_b, 'out1', evaporator, 'in1')
evaporator_sink_b = connection(evaporator, 'out1', preheater, 'in1')
preheater_sink = connection(preheater, 'out1', sink_b, 'in1')
nw.add_conns(evaporator_brine_in, evaporator_sink_b, preheater_sink)
# cooling air cycle
ca_in = connection(source_ca, 'out1', condenser, 'in2')
ca_out = connection(condenser, 'out2', sink_ca, 'in1')
nw.add_conns(ca_in, ca_out)

# parametrization of components
evaporator.set_attr(pr1=0.93181818, pr2=0.970588)
preheater.set_attr(pr1=0.949494, pr2=0.955752)
turbine.set_attr(pr=0.098148148, eta_s=0.85, design=['eta_s', 'pr'])
pump.set_attr(eta_s=0.9)
ihe.set_attr(pr1=0.849056603, pr2=0.957627118)
condenser.set_attr(pr1=0.8889, pr2=1)
ihe.set_attr(ttd_u=16.7)

# busses
# characteristic function for generator efficiency
# x = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
# y = np.array([0, 0.88, 0.89, 0.90, 0.91, 0.976, 0.91])
# gen = char_line(x=x, y=y)
# motor of pump has a constant efficiency
# power = bus('total output power')
# power.add_comps({'c': turbine, 'p': 'P', 'char': gen})
# nw.add_busses(power)

# parametrization of connections
evaporator_turbine.set_attr(p=p_before_turbine, state='g', fluid={'water': 0, 'Isopentane': 1, 'Air': 0})

evaporator_brine_in.set_attr(T=T_brine_in, m=mass_flow_rate_brine, fluid={'water': 1, 'Isopentane': 0, 'Air':0})
evaporator_sink_b.set_attr(T=T_brine_in-28)
preheater_sink.set_attr(T=T_reinjection)

# air cooling connections
ca_in.set_attr(T=T_air, p=p_air, fluid={'water': 0, 'Isopentane': 0, 'Air': 1})
ca_out.set_attr(T=T_air + 15)

# solving
mode = 'design'
save_path = 'power_generation_with_low_T'
# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
nw.save(save_path)