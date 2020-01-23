from tespy.connections import connection, ref
from tespy.networks import network
from tespy.components import (condenser, heat_exchanger, pump, turbine, source,
                              sink, cycle_closer, splitter, merge, valve)
from CoolProp.CoolProp import PropsSI
import numpy as np
from tespy.tools import logger
import logging
mypath = logger.define_logging(
log_path=True, log_version=True, timed_rotating={'backupCount': 4},
screen_level=logging.WARNING, screen_datefmt = "no_date")
# input parameters
fluids = ['water', 'Isopentane', 'Air']

p_before_turbines = 11.49

mass_flow_rate_air = 6142
T_air = -4.7
p_air = 0.61

mass_flow_rate_brine = 3.4199e2
T_brine_in = 146.6
p_brine_in = 9.4
T_brine_out = 69.1

# basic network
nw = network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

# main components
cond1 = condenser('condenser 1')
cond2 = condenser('condenser 2')
ihe1 = heat_exchanger('internal heat exchanger 1')
ihe2 = heat_exchanger('internal heat exchanger 2')
pu1 = pump('feeding pump 1')
pu2 = pump('feeding pump 2')
turb1 = turbine('turb 1')
turb2 = turbine('turb 2')
val1 = valve('control valve 1')
val2 = valve('control valve 2')

spl = splitter('main cycle splitter')
mer = merge('main cycle merge')

preh = heat_exchanger('preheater')
evap = heat_exchanger('evaporator')

# cooling air
source_ca_1 = source('cooling air source 1')
sink_ca_1 = sink('cooling air sink 1')
source_ca_2 = source('cooling air source 2')
sink_ca_2 = sink('cooling air sink 2')

#brine
source_b = source('brine source')
sink_b = sink('brine sink')

# working fluid
close_cycle = cycle_closer('cycle closer before turbines')

# connections
wf_in_spl = connection(close_cycle, 'out1', spl, 'in1')
nw.add_conns(wf_in_spl)

# main cycle part 1
spl_turb_1 = connection(spl, 'out1', turb1, 'in1')
turb_ihe_1 = connection(turb1, 'out1', ihe1, 'in1')
ihe_cond_1 = connection(ihe1, 'out1', cond1, 'in1')
cond_pu_1 = connection(cond1, 'out1', pu1, 'in1')
pu_ihe_1 = connection(pu1, 'out1', ihe1, 'in2')
ihe_val_1 = connection(ihe1, 'out2', val1, 'in1')
val_mer_1 = connection(val1, 'out1', mer, 'in1')
nw.add_conns(spl_turb_1, turb_ihe_1, ihe_cond_1, cond_pu_1,
             pu_ihe_1, ihe_val_1, val_mer_1)

# main cycle part 2
spl_turb_2 = connection(spl, 'out2', turb2, 'in1')
turb_ihe_2 = connection(turb2, 'out1', ihe2, 'in1')
ihe_cond_2 = connection(ihe2, 'out1', cond2, 'in1')
cond_pu_2 = connection(cond2, 'out1', pu2, 'in1')
pu_ihe_2 = connection(pu2, 'out1', ihe2, 'in2')
ihe_val_2 = connection(ihe2, 'out2', val2, 'in1')
val_mer_2 = connection(val2, 'out1', mer, 'in2')
nw.add_conns(spl_turb_2, turb_ihe_2, ihe_cond_2, cond_pu_2,
             pu_ihe_2, ihe_val_2, val_mer_2)

# preheater and evaporator
mer_preh = connection(mer, 'out1', preh, 'in2')
preh_evap = connection(preh, 'out2', evap, 'in2')
evap_wf_out = connection(evap, 'out2', close_cycle, 'in1')
nw.add_conns(mer_preh, preh_evap, evap_wf_out)

# geo-fluid
brine_in = connection(source_b, 'out1', evap, 'in1')
evap_preh = connection(evap, 'out1', preh, 'in1')
brine_out = connection(preh, 'out1', sink_b, 'in1')
nw.add_conns(brine_in, evap_preh, brine_out)

# cooling air cycle 1
ca_in_1 = connection(source_ca_1, 'out1', cond1, 'in2')
ca_out_1 = connection(cond1, 'out2', sink_ca_1, 'in1')
nw.add_conns(ca_in_1, ca_out_1)

# cooling air cycle 2
ca_in_2 = connection(source_ca_2, 'out1', cond2, 'in2')
ca_out_2 = connection(cond2, 'out2', sink_ca_2, 'in1')
nw.add_conns(ca_in_2, ca_out_2)

# parametrization of components main cycle part 1
turb1.set_attr()
ihe1.set_attr(pr1=0.849056603, pr2=0.957627118)
cond1.set_attr(pr1=0.8889, pr2=1)
pu1.set_attr(eta_s=0.9, pr=14)

# parametrization of components main cycle part 2
turb2.set_attr()
ihe2.set_attr(pr1=0.849056603, pr2=0.957627118)
cond2.set_attr(pr1=0.8889, pr2=1)
pu2.set_attr(eta_s=0.9, pr=14)

# parametrization of preheater and evaporator
preh.set_attr(pr1=1, pr2=0.955752212)
evap.set_attr(pr1=1, pr2=1, ttd_l=10)

# parametrization of connections main cycle part 1
spl_turb_1.set_attr()
turb_ihe_1.set_attr(p=1.31, T=76)
ihe_cond_1.set_attr()
cond_pu_1.set_attr()
pu_ihe_1.set_attr()
ihe_val_1.set_attr()
val_mer_1.set_attr()

# parametrization of connections main cycle part 2
spl_turb_2.set_attr()
turb_ihe_2.set_attr(p=1.29, T=71.6)
ihe_cond_2.set_attr()
cond_pu_2.set_attr()
pu_ihe_2.set_attr()
ihe_val_2.set_attr()
val_mer_2.set_attr()

# parametrization of connections (preheater and evaporator)
wf_in_spl.set_attr(Td_bp=2.3, h0=500, p=p_before_turbines, fluid={'water': 0, 'Isopentane': 1, 'Air': 0})
spl_turb_2.set_attr(m=ref(wf_in_spl, 0.5, 0))
evap_preh.set_attr(T=80)

# air cooling connections
ca_in_1.set_attr(T=T_air, p=p_air, m=mass_flow_rate_air / 2, fluid={'water': 0, 'Isopentane': 0, 'Air': 1})
ca_out_1.set_attr(T=T_air + 15)

ca_in_2.set_attr(T=T_air, p=p_air, m=mass_flow_rate_air / 2, fluid={'water': 0, 'Isopentane': 0, 'Air': 1})
ca_out_2.set_attr(T=T_air + 15)

# brine connections
brine_in.set_attr(T=T_brine_in, p=p_brine_in, state='l', fluid={'water': 1, 'Isopentane': 0, 'Air': 0})
brine_out.set_attr(T=T_brine_out)
# solving
mode = 'design'
file = 'yangyi'
# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
nw.save(file)
