from tespy.connections import connection, ref
from tespy.networks import network
from tespy.components import condenser, heat_exchanger, pump, turbine, source, sink, cycle_closer
from CoolProp.CoolProp import PropsSI
import numpy as np
from tespy.tools import logger
import logging
mypath = logger.define_logging(
log_path=True, log_version=True, timed_rotating={'backupCount': 4},
screen_level=logging.WARNING, screen_datefmt = "no_date")
# input parameters
fluids = ['water', 'Isopentane', 'Air']

mass_flow_rate_brine = 3.4199e2
T_brine_in = 146.6
p_brine_in = 9.4
T_brine_out = 69.1

mass_flow_rate_air = 6142
T_air = -4.7
p_air = 0.61


# calculation secondary variables
p_before_turbine = PropsSI('P', 'T', T_brine_in + 273.15 - 26, 'Q', 1, 'Isopentane')/1e5
# basic network
nw = network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

# main components
condenser = condenser('condenser')
ihe = heat_exchanger('internal heat exchanger')
pump = pump('feeding pump')
turbine = turbine('turbine')
preh = heat_exchanger('preheater')
evap = heat_exchanger('evaporator')

# cooling air
source_ca = source('cooling air source')
sink_ca = sink('cooling air sink')

#brine
source_b = source('brine source')
sink_b = sink('brine sink')

# working fluid
close_cycle = cycle_closer('cycle closer before turbine')

# connections
# main cycle
turbine_wf_in = connection(close_cycle, 'out1', turbine, 'in1')
turbine_ihe = connection(turbine, 'out1', ihe, 'in1')
ihe_condenser = connection(ihe, 'out1', condenser, 'in1')
condenser_pump = connection(condenser, 'out1', pump, 'in1')
pump_ihe = connection(pump, 'out1', ihe, 'in2')
ihe_preh = connection(ihe, 'out2', preh, 'in2')
preh_evap = connection(preh, 'out2', evap, 'in2')
evap_wf_out = connection(evap, 'out2', close_cycle, 'in1')
nw.add_conns(turbine_wf_in, turbine_ihe, ihe_condenser, condenser_pump, pump_ihe, ihe_preh, preh_evap, evap_wf_out)

# geo-fluid
brine_in = connection(source_b, 'out1', evap, 'in1')
evap_preh = connection(evap, 'out1', preh, 'in1')
brine_out = connection(preh, 'out1', sink_b, 'in1')
nw.add_conns(brine_in, evap_preh, brine_out)

# cooling air
ca_in = connection(source_ca, 'out1', condenser, 'in2')
ca_out = connection(condenser, 'out2', sink_ca, 'in1')
nw.add_conns(ca_in, ca_out)

# busses
# characteristic function for generator efficiency
#x = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
#y = np.array([0, 0.88, 0.89, 0.90, 0.91, 0.976, 0.91])
#gen = cmp_char.characteristics(x=x, y=y)
## motor of pump has a constant efficiency
#power = bus('total output power')
#power.add_comps({'c': turbine, 'p': 'P', 'char': gen},
#                {'c': pump, 'char': 1 / 0.85})
#nw.add_busses(power)
# parametrization of components
#v = np.array([0.1, 0.4])
#dp = np.array([11, 11]) * 1e5
#char = hlp.dc_cc(x=v, y=dp, is_set=True)
#pump.set_attr(eta_s=0.85, flow_char=char, design=['eta_s_char'])

# parametrization of components
pump.set_attr(eta_s=0.9)
ihe.set_attr(pr1=0.849056603, pr2=0.957627118)
condenser.set_attr(pr1=0.8889, pr2=1)
turbine.set_attr(pr=0.098148148, eta_s=0.85, design=['eta_s', 'pr'])
preh.set_attr(pr1=1, pr2=0.955752212)
evap.set_attr(pr1=1, pr2=1, ttd_l=5)

# parametrization of connections
turbine_wf_in.set_attr(Td_bp=2.3, p=p_before_turbine, h0=500, fluid={'water': 0, 'Isopentane': 1, 'Air': 0})
#condenser_pump.set_attr(Td_bp=0)

# air cooling connections
ca_in.set_attr(T=T_air, p=p_air, m=mass_flow_rate_air, fluid={'water': 0, 'Isopentane': 0, 'Air': 1})
ca_out.set_attr(T=T_air + 15)

# brine connections
brine_in.set_attr(T=T_brine_in, p=p_brine_in, m=mass_flow_rate_brine, state='l', fluid={'water': 1, 'Isopentane': 0, 'Air': 0})
brine_out.set_attr(T=T_brine_out)
# solving
mode = 'design'
file = 'yangyi'
# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
#nw.save(file)
