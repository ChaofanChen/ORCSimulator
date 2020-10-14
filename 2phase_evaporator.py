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
mass_flow_rate_brine = 180
mass_flow_rate_steam = 5
T_brine_in = 140
T_reinjection = 70
# calculation secondary variables
p_steam_in = PropsSI('P', 'T', T_brine_in+273.15, 'Q', 1, 'water')/1e5
# main components
evaporator = orc_evaporator('orc_evaporator')
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
wf_in_evaporator = connection(source_wf, 'out1', evaporator, 'in3')
evaporator_wf_out = connection(evaporator, 'out3', sink_wf, 'in1')
nw.add_conns(wf_in_evaporator, evaporator_wf_out)
# geo-steam cycle
evaporator_steam_in = connection(source_s, 'out1', evaporator, 'in1')
evaporator_sink_s = connection(evaporator, 'out1', sink_s, 'in1')
nw.add_conns(evaporator_steam_in, evaporator_sink_s)
# geo-brine cycle
evaporator_brine_in = connection(source_b, 'out1', evaporator, 'in2')
evaporator_sink_b = connection(evaporator, 'out2', sink_b, 'in1')
nw.add_conns(evaporator_brine_in, evaporator_sink_b)
# parametrization of components
evaporator.set_attr(pr1=0.90, pr2=0.970588, pr3=1)
evaporator.set_attr(ttd_u=20.8, design=['ttd_u'])
# parametrization of connections
evaporator_steam_in.set_attr(T=T_brine_in, m=mass_flow_rate_steam, p=p_steam_in, state='g', fluid={'water': 1, 'Isopentane': 0, 'Air':0})
evaporator_brine_in.set_attr(T=T_brine_in, p=p_steam_in, m=mass_flow_rate_brine, state='l', fluid={'water': 1, 'Isopentane': 0, 'Air':0})
wf_in_evaporator.set_attr(T=100.3, m=43, fluid={'water': 0, 'Isopentane': 1, 'Air':0})
# solving
mode = 'design'
nw.solve(mode=mode)
nw.print_results()
