from tespy import con, nwk, cmp, hlp, cmp_char
from CoolProp.CoolProp import PropsSI
import numpy as np
from tespy.tools import logger
import logging
mypath = logger.define_logging(
log_path=True, log_version=True, timed_rotating={'backupCount': 4},
screen_level=logging.WARNING, screen_datefmt = "no_date")
# input parameters
fluids = ['water', 'Isopentane', 'N2', 'Ar', 'O2']

p_after_pump = 11.8

mass_flow_rate_air = 6142
T_air = -4.7
p_air = 0.61

mass_flow_rate_brine = 3.4199e2
T_brine_in = 146.6
p_brine_in = 9.4
T_brine_out = 69.1

# calculation
T_before_turbine = PropsSI('T', 'P', p_after_pump*0.957627118*0.955752212*1e5, 'Q', 1, 'Isopentane')-273.15+2.3
# basic network
nw = nwk.network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
# main components
condenser = cmp.condenser('condenser')
ihe = cmp.heat_exchanger('internal heat exchanger')
pump = cmp.pump('feeding pump')
turbine = cmp.turbine('turbine')
p_and_e = cmp.heat_exchanger('preheater and evaporator')
# cooling air
source_ca = cmp.source('cooling air source')
sink_ca = cmp.sink('cooling air sink')
#brine
source_b = cmp.source('brine source')
sink_b = cmp.sink('brine sink')
# working fluid
source_wf_1 = cmp.source('working fluid source before turbine')
sink_wf_1 = cmp.sink('working fluid sink from before turbine')
source_wf_2 = cmp.source('working fluid source from ihe')
sink_wf_2 = cmp.sink('working fluid sink from ihe')
# connections
# main cycle
p_and_e_wf_in = con.connection(source_wf_2, 'out1', p_and_e, 'in2')
p_and_e_wf_out = con.connection(p_and_e, 'out2', sink_wf_2, 'in1')
turbine_wf_in = con.connection(source_wf_1, 'out1', turbine, 'in1')
turbine_ihe = con.connection(turbine, 'out1', ihe, 'in1')
ihe_condenser = con.connection(ihe, 'out1', condenser, 'in1')
condenser_pump = con.connection(condenser, 'out1', pump, 'in1')
pump_ihe = con.connection(pump, 'out1', ihe, 'in2')
ihe_wf_out = con.connection(ihe, 'out2', sink_wf_1, 'in1')
nw.add_conns(turbine_wf_in, turbine_ihe, ihe_condenser, condenser_pump, pump_ihe, ihe_wf_out)
# geo-fluid
brine_in = con.connection(source_b, 'out1', p_and_e, 'in1')
brine_out = con.connection(p_and_e, 'out1', sink_b, 'in1')
nw.add_conns(p_and_e_wf_in, p_and_e_wf_out, brine_in, brine_out)
# cooling air
ca_in = con.connection(source_ca, 'out1', condenser, 'in2')
ca_out = con.connection(condenser, 'out2', sink_ca, 'in1')
nw.add_conns(ca_in, ca_out)
# busses
# characteristic function for generator efficiency
#x = np.array([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
#y = np.array([0, 0.88, 0.89, 0.90, 0.91, 0.976, 0.91])
#gen = cmp_char.characteristics(x=x, y=y)
## motor of pump has a constant efficiency
#power = con.bus('total output power')
#power.add_comps({'c': turbine, 'p': 'P', 'char': gen},
#                {'c': pump, 'char': 1 / 0.85})
#nw.add_busses(power)
# parametrization of components
#v = np.array([0.1, 0.4])
#dp = np.array([11, 11]) * 1e5
#char = hlp.dc_cc(x=v, y=dp, is_set=True)
#pump.set_attr(eta_s=0.85, flow_char=char, design=['eta_s_char'])
pump.set_attr(pr=14.75, eta_s=0.9)
ihe.set_attr(pr1=0.849056603, pr2=0.957627118)
condenser.set_attr(pr1=0.8889, pr2=1)
turbine.set_attr(pr=0.098148148, eta_s=0.85, design=['eta_s', 'pr'])
p_and_e.set_attr(pr1=1, pr2=0.955752212)
# parametrization of connections
turbine_wf_in.set_attr(T=T_before_turbine, state='g', 
                       fluid={'water': 0, 'Isopentane': 1, 'N2': 0, 'Ar': 0, 'O2': 0})
p_and_e_wf_in.set_attr(T=con.ref(ihe_wf_out, 1, 0), p=con.ref(ihe_wf_out, 1, 0), 
                       m=con.ref(ihe_wf_out, 1, 0), fluid={'water': 0, 'Isopentane': 1, 'N2': 0, 'Ar': 0, 'O2': 0})
p_and_e_wf_out.set_attr(T=con.ref(turbine_wf_in, 1, 0), state='g')
pump_ihe.set_attr(p=p_after_pump)
# air cooling connections
ca_in.set_attr(T=T_air, p=p_air, m=mass_flow_rate_air,
               fluid={'water': 0, 'Isopentane': 0, 'N2': 0.7812, 'Ar': 0.0092, 'O2': 0.2096})
ca_out.set_attr(T=T_air+15)
# brine connections
brine_in.set_attr(T=T_brine_in, p=p_brine_in, m=mass_flow_rate_brine, state='l', 
                  fluid={'water': 1, 'Isopentane': 0, 'N2': 0, 'Ar': 0, 'O2': 0})
brine_out.set_attr(T=T_brine_out)
# solving
mode = 'design'
file = 'yangyi'
# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
nw.save(file, structure=True)