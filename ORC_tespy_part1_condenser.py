from tespy.connections import connection
from tespy.networks import network
from tespy.components import condenser, source, sink
from CoolProp.CoolProp import PropsSI
import numpy as np
from tespy.tools import logger
import logging
mypath = logger.define_logging(
    log_path=True, log_version=True, timed_rotating={'backupCount': 4},
    screen_level=logging.WARNING, screen_datefmt = "no_date")
fluids = ['Isopentane', 'air']
T_air = -4.7
p_air = 0.61
nw = network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg',
            p_range=[0.01, 150], T_range=[5, 800], h_range=[-10, 5000])
# components
# main components
condenser = condenser('condenser')
# cooling air
source_ca = source('cooling air source')
sink_ca = sink('cooling air sink')
# working fluid
source_wf = source('working fluid source')
sink_wf = sink('working fluid sink')
# connections
# main cycle
# geo-fluid
ihe_condenser = connection(source_wf, 'out1', condenser, 'in1')
condenser_pump = connection(condenser, 'out1', sink_wf, 'in1')
nw.add_conns(ihe_condenser, condenser_pump)
# cooling air
ca_in = connection(source_ca, 'out1', condenser, 'in2')
ca_out = connection(condenser, 'out2', sink_ca, 'in1')
nw.add_conns(ca_in, ca_out)
# parametrization of components
condenser.set_attr(pr1=0.8889, pr2=1)
#condenser.set_attr(Q=-9.365e7)
# parametrization of connections
ihe_condenser.set_attr(T=42.4, p=0.9, m=243.72, fluid={'Isopentane': 1, 'air': 0})

ca_in.set_attr(T=T_air, p=p_air, fluid={'Isopentane': 0, 'air': 1})
ca_out.set_attr(T=T_air+15)
# solving
mode = 'design'
file = 'yangyi_condenser'
# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
#nw.save(file)

