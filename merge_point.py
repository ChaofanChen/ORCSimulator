from tespy.connections import connection
from tespy.networks import network
from tespy.components import evaporator, heat_exchanger, pump, turbine, source, sink, cycle_closer, splitter, merge
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
pump = pump('condensate pump')
merge = merge('geo-fluid merge point')
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

scr_s_merge = connection(source_s, 'out1', merge, 'in1')

scr_b_merge = connection(source_b, 'out1', merge, 'in2')

merge_sink = connection(merge, 'out1', sink_b, 'in1')

nw.add_conns(scr_s_merge, scr_b_merge, merge_sink)
# parametrization of connections
scr_s_merge.set_attr(T=40, p=10, m=20.4, fluid={'water': 1, 'Isopentane': 0})
scr_b_merge.set_attr(T=40, m=190, fluid={'water': 1, 'Isopentane': 0})
#merge_sink.set_attr(fluid={'water': 1, 'Isopentane': 0})
# solving
mode = 'design'
file = 'yangyi_evaporator_new'
# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
#nw.save(file)
