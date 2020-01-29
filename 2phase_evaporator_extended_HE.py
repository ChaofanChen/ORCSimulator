from tespy.connections import connection
from tespy.networks import network
from tespy.components import evaporator, heat_exchanger, pump, turbine, source, sink, cycle_closer
from CoolProp.CoolProp import PropsSI
import numpy as np
from tespy.tools import logger
import logging
mypath = logger.define_logging(
log_path=True, log_version=True, timed_rotating={'backupCount': 4},
screen_level=logging.WARNING, screen_datefmt = "no_date")

fluids = ['water', 'Isopentane']

nw = network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
# components
# main components
evaporator = evaporator('evaporator')
# working fluid
source_wf = source('working fluid source')
sink_wf = sink('working fluid sink')
#brine
source_s = source('steam source')
source_b = source('brine source')
sink_s = sink('steam sink')
sink_b = sink('brine sink')
# connections
# main cycle
evaporator_wf_in = connection(source_wf, 'out1', evaporator, 'in3')
evaporator_wf_out = connection(evaporator, 'out3', sink_wf, 'in1')

evaporator_steam_in = connection(source_s, 'out1', evaporator, 'in1')
evaporator_sink_s = connection(evaporator, 'out1', sink_s, 'in1')

evaporator_brine_in = connection(source_b, 'out1', evaporator, 'in2')
evaporator_sink_b = connection(evaporator, 'out2', sink_b, 'in1')

nw.add_conns(evaporator_wf_in, evaporator_wf_out)
nw.add_conns(evaporator_steam_in, evaporator_sink_s)
nw.add_conns(evaporator_brine_in, evaporator_sink_b)
# parametrization of components
evaporator.set_attr(pr1=0.93181818, pr2=0.970588, pr3=1)

# parametrization of connections
evaporator_wf_in.set_attr(T=111.6, p=10.8, m=243.72, fluid={'water': 0, 'Isopentane': 1})
evaporator_wf_out.set_attr(T=119.8)

evaporator_steam_in.set_attr(T=146.6, p=4.4, m=21, fluid={'water': 1, 'Isopentane': 0})
evaporator_sink_s.set_attr(T=132.5)

evaporator_brine_in.set_attr(T=146.6, p=10.2, fluid={'water': 1, 'Isopentane': 0})
evaporator_sink_b.set_attr(T=118.6)
# solving
mode = 'design'
file = 'yangyi_evaporator_new'
# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
#nw.save(file)
