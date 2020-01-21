from tespy import con, nwk, cmp, hlp, cmp_char
#from matplotlib import pyplot as plt
import numpy as np
from tespy.tools import logger
import logging
mypath = logger.define_logging(
log_path=True, log_version=True, timed_rotating={'backupCount': 4},
screen_level=logging.WARNING, screen_datefmt = "no_date")
#

fluids = ['water', 'Isopentane']

nw = nwk.network(fluids=fluids)
nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')
nw.set_printoptions(print_level='info')
# components
# main components
preheater = cmp.heat_exchanger('preheater')
# working fluid
source_wf = cmp.source('working fluid source')
sink_wf = cmp.sink('working fluid sink')
#brine
source_b = cmp.source('brine source')
sink_b = cmp.sink('brine sink')
# connections
# main cycle
# geo-fluid
preheater_wf_in = con.connection(source_wf, 'out1', preheater, 'in2', m=243.72)
preheater_wf_out = con.connection(preheater, 'out2', sink_wf, 'in1')
brine_in = con.connection(source_b, 'out1', preheater, 'in1')
brine_out = con.connection(preheater, 'out1', sink_b, 'in1')
nw.add_conns(preheater_wf_in, preheater_wf_out, brine_in, brine_out)
# parametrization of components
preheater.set_attr(pr1=0.949494949494, pr2=0.955752212)
#p_and_e.set_attr(Q=-9.365e7)
# parametrization of connections
preheater_wf_in.set_attr(T=39.7, p=11.3, fluid={'water': 0, 'Isopentane': 1})
preheater_wf_out.set_attr(T=111.6)

brine_in.set_attr(T=119.9, fluid={'water': 1, 'Isopentane': 0})
brine_out.set_attr(T=69.1, p=9.4)
# solving
mode = 'design'
file = 'yangyi_preheater'
# solve the network, print the results to prompt and save
nw.solve(mode=mode)
nw.print_results()
nw.save(file)
