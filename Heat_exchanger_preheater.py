from tespy.connections import connection
from tespy.networks import network
from tespy.components import heat_exchanger, source, sink
import numpy as np
from tespy.tools import logger
import logging
mypath = logger.define_logging(
log_path=True, log_version=True, timed_rotating={'backupCount': 4},
screen_level=logging.WARNING, screen_datefmt = "no_date")
#

fluids = ['water', 'Isopentane']

nw = network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
# components
# main components
preheater = heat_exchanger('preheater')
# working fluid
source_wf = source('working fluid source')
sink_wf = sink('working fluid sink')
#brine
source_b = source('brine source')
sink_b = sink('brine sink')
# connections
# main cycle
# geo-fluid
preheater_wf_in = connection(source_wf, 'out1', preheater, 'in2')
preheater_wf_out = connection(preheater, 'out2', sink_wf, 'in1')
brine_in = connection(source_b, 'out1', preheater, 'in1')
brine_out = connection(preheater, 'out1', sink_b, 'in1')

nw.add_conns(preheater_wf_in, preheater_wf_out)
nw.add_conns(brine_in, brine_out)
#nw.add_conns(preheater_wf_in, preheater_wf_out, brine_in, brine_out)
# parametrization of components
preheater.set_attr(pr1=0.949494949494, pr2=0.955752212, ttd_l=29.4)
#p_and_e.set_attr(Q=-9.365e7)
# parametrization of connections
preheater_wf_in.set_attr(T=39.7, p=11.3, m=243.72, fluid={'water': 0, 'Isopentane': 1})
preheater_wf_out.set_attr(T=111.6)

brine_in.set_attr(T=119.9, p=9.9, fluid={'water': 1, 'Isopentane': 0})
# brine_out.set_attr(T=69.1)
# solving
mode = 'design'
file = 'yangyi_preheater'
# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
#nw.save(file)
