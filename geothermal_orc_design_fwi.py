from tespy.connections import connection, bus, ref
from tespy.tools import char_line
from tespy.networks import network
from tespy.components import (
    heat_exchanger, pump, turbine, source, sink, cycle_closer, splitter,
    merge, condenser, drum, valve, heat_exchanger_simple
)
from tespy.tools import logger

from fluprodia import FluidPropertyDiagram
from CoolProp.CoolProp import PropsSI as PSI

import pandas as pd
import numpy as np
import logging


logger.define_logging(screen_level=logging.ERROR)


class PowerPlant():

    def __init__(self, working_fluid):

        self.working_fluid = working_fluid
        fluids = ['water', self.working_fluid, 'air']
        self.nw = network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # geo parameters

        geo_mass_flow = 210
        geo_steam_share = 0.1
        T_brine_in = 144.8
        T_reinjection = 70.8

        # ambient parameters

        T_amb = 0
        p_amb = 1

        # main components

        geo_steam = source('steam source')
        geo_brine = source('brine source')
        geo_reinjection = sink('reinjection well')

        air_in = source('ambient air source')
        air_out = sink('ambient air sink')
        air_cond = condenser('main condenser')

        orc_cc = cycle_closer('orc cycle closer')

        evap_steam = condenser('steam evaporator')
        evap_splitter = splitter('splitter evaporation')
        evap_merge = merge('merge evaporation')
        evap_steam = condenser('steam evaporator')
        geo_steam_pump = pump('geosteam condensate pump')
        evap_brine = heat_exchanger('brine evaporator')
        dr = drum('drum')

        eco = heat_exchanger('economiser')
        feed_water_pump = pump('feed water pump')
        geo_merge = merge('brine merge')

        tur = turbine('turbine')

        ls_valve = valve('live steam valve')

        ihe = heat_exchanger('internal heat exchanger')

        # busses
        power_bus = bus('power output')
        power_bus.add_comps(
            {'c': tur, 'char': -1},
            {'c': feed_water_pump, 'char': -1}, {'c': geo_steam_pump, 'char': -1}
        )

        geothermal_bus = bus('thermal input')
        geothermal_bus.add_comps(
            {'c': eco, 'char': -1}, {'c': evap_brine, 'char': -1},
            {'c': evap_steam, 'char': -1}
        )

        self.nw.add_busses(power_bus, geothermal_bus)

        # turbine to condenser
        ls_in = connection(orc_cc, 'out1', ls_valve, 'in1')
        lsv_tur = connection(ls_valve, 'out1', tur, 'in1')
        tur_ihe = connection(tur, 'out1', ihe, 'in1')
        ihe_cond = connection(ihe, 'out1', air_cond, 'in1')
        self.nw.add_conns(ls_in, lsv_tur, tur_ihe, ihe_cond)

        # condenser to steam generator
        cond_fwp = connection(air_cond, 'out1', feed_water_pump, 'in1')
        fwp_ihe = connection(feed_water_pump, 'out1', ihe, 'in2')
        self.nw.add_conns(cond_fwp, fwp_ihe)

        # steam generator
        ihe_eco = connection(ihe, 'out2', eco, 'in2')
        eco_dr = connection(eco, 'out2', dr, 'in1')
        dr_esp = connection(dr, 'out1', evap_splitter, 'in1')
        esp_eb = connection(evap_splitter, 'out1', evap_brine, 'in2')
        esp_es = connection(evap_splitter, 'out2', evap_steam, 'in2')
        eb_em = connection(evap_brine, 'out2', evap_merge, 'in1')
        es_em = connection(evap_steam, 'out2', evap_merge, 'in2')
        em_dr = connection(evap_merge, 'out1', dr, 'in2')
        ls_out = connection(dr, 'out2', orc_cc, 'in1')
        self.nw.add_conns(ihe_eco, eco_dr, dr_esp, esp_eb, esp_es, eb_em, es_em, em_dr, ls_out)

        # air cold side
        air_cold = connection(air_in, 'out1', air_cond, 'in2')
        air_hot = connection(air_cond, 'out2', air_out, 'in1')
        self.nw.add_conns(air_cold, air_hot)

        # geo source
        gs_es = connection(geo_steam, 'out1', evap_steam, 'in1', label='geosteam')
        es_gsp = connection(evap_steam, 'out1', geo_steam_pump, 'in1')
        gsp_gm = connection(geo_steam_pump, 'out1', geo_merge, 'in1')
        gb_eb = connection(geo_brine, 'out1', evap_brine, 'in1', label='geobrine')
        eb_gm = connection(evap_brine, 'out1', geo_merge, 'in2')
        self.nw.add_conns(gs_es, es_gsp, gsp_gm, gb_eb, eb_gm)

        gm_eco = connection(geo_merge, 'out1', eco, 'in1')
        eco_gr = connection(eco, 'out1', geo_reinjection, 'in1', label='reinjection')
        self.nw.add_conns(gm_eco, eco_gr)

        # fluid settings
        ihe_eco.set_attr(fluid={self.working_fluid: 1, 'air': 0, 'water': 0})
        air_cold.set_attr(fluid={self.working_fluid: 0, 'air': 1, 'water': 0})
        gs_es.set_attr(fluid={self.working_fluid: 0, 'air': 0, 'water': 1})
        gb_eb.set_attr(fluid={self.working_fluid: 0, 'air': 0, 'water': 1})

        # connection parameters
        ls_stable_p0 = PSI('P', 'T', T_brine_in + 273.15, 'Q', 1, self.working_fluid) / 1e5
        lsv_tur.set_attr(p0=ls_stable_p0)
        ws_stable_h0 = (
            PSI('H', 'T', T_amb + 273.15, 'Q', 1, self.working_fluid) + 0.5 * (
                PSI('H', 'T', T_brine_in + 273.15, 'Q', 1, self.working_fluid) -
                PSI('H', 'T', T_amb + 273.15, 'Q', 1, self.working_fluid)
            )
        ) / 1e3
        tur_ihe.set_attr(h=ws_stable_h0)
        ihe_cond.set_attr(Td_bp=2, design=['Td_bp'], p0=PSI('P', 'T', T_amb + 273.15, 'Q', 1, self.working_fluid) / 1e5)
        fwp_ihe.set_attr(h=ref(cond_fwp, 1, 1e3))

        # steam generator
        gs_es.set_attr(m=geo_mass_flow * geo_steam_share, T=T_brine_in, x=1, p0=5)
        gb_eb.set_attr(m=geo_mass_flow * (1 - geo_steam_share), T=T_brine_in, x=0)

        em_dr.set_attr()
        eb_em.set_attr(x=0.5)
        es_em.set_attr(x=0.5, design=['x'])
        eb_gm.set_attr(T=T_brine_in - 20)

        eco_dr.set_attr(Td_bp=-2)

        # main condenser
        air_cold.set_attr(p=p_amb, T=T_amb)
        air_hot.set_attr(T=15)

        # component parameters
        # turbines
        tur.set_attr(design=['eta_s'], offdesign=['cone', 'eta_s_char'])
        ls_valve.set_attr(pr=1, design=['pr'])
        # condensing
        ihe.set_attr(pr1=1, pr2=1, offdesign=['kA_char'])
        air_cond.set_attr(pr1=1, pr2=1, ttd_u=10, design=['ttd_u'], offdesign=['kA_char'])
        feed_water_pump.set_attr(design=['eta_s'], offdesign=['eta_s_char'])

        # steam generator
        evap_steam.set_attr(pr1=0.99, offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        evap_brine.set_attr(pr1=1, ttd_l=10, offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        eco.set_attr(pr1=1, pr2=1)
        geo_steam_pump.set_attr(eta_s=0.75, design=['eta_s'], offdesign=['eta_s_char'])

        self.nw.set_attr(iterinfo=False)
        self.nw.solve('design')
        self.nw.print_results()
        tur.set_attr(eta_s=0.9)
        feed_water_pump.set_attr(eta_s=0.75)
        tur_ihe.set_attr(h=None)
        fwp_ihe.set_attr(h=None)
        eb_gm.set_attr(T=None)

    def calculate_efficiency(self, geo_mass_flow, geo_steam_fraction, T_reinjection):

        self.nw.connections['geosteam'].set_attr(m=geo_mass_flow * geo_steam_fraction)
        self.nw.connections['geobrine'].set_attr(m=geo_mass_flow * (1 - geo_steam_fraction))
        self.nw.connections['reinjection'].set_attr(T=T_reinjection)
        self.nw.solve('design')

        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            return np.nan
        else:
            return self.nw.busses['power output'].P.val / self.nw.busses['thermal input'].P.val

    def print_result(self):
        eta_th = self.nw.busses['power output'].P.val / self.nw.busses['thermal input'].P.val
        power = self.nw.busses['power output'].P.val
        print('Power output: {} W'.format(round(power, 0)))
        print('Thermal efficiency: {} %'.format(round(eta_th * 100, 2)))

    def plot_process(self, fn='somename'):

        result_dict = {
            prop: [
                conn.get_attr(prop).val for conn in self.nw.conns.index
                if conn.fluid.val[self.working_fluid] == 1
            ] for prop in ['p', 'h', 's', 'T']
        }

        self.diagram.set_limits(
            x_min=min(result_dict['h']) - 50,
            x_max=max(result_dict['h']) + 50,
            y_min=min(result_dict['p']) / 2,
            y_max=max(result_dict['p']) * 10
        )
        self.diagram.draw_isolines('logph')
        self.diagram.ax.scatter(result_dict['h'], result_dict['p'], zorder=100)
        self.diagram.save(fn + '_logph.pdf')

        self.diagram.set_limits(
            x_min=min(result_dict['s']) - 50,
            x_max=max(result_dict['s']) + 50,
            y_min=min(result_dict['T']) - 25,
            y_max=PSI('T_critical', self.working_fluid) + 25 - 273.15
        )
        self.diagram.draw_isolines('Ts')
        self.diagram.ax.scatter(result_dict['s'], result_dict['T'], zorder=100)
        self.diagram.save(fn + '_Ts.pdf')

    def generate_diagram(self):

        self.diagram = FluidPropertyDiagram(self.working_fluid)
        self.diagram.set_unit_system(T='°C', p='bar', h='kJ/kg')
        iso_T = np.arange(0, 201, 25)
        self.diagram.set_isolines(T=iso_T)
        self.diagram.calc_isolines()


fluids = ['R600', 'R245fa', 'R245CA', 'Toluene', 'Isopentane', 'n-Pentane', 'R123']
for fluid in fluids:
    try:
        # for some testing
        sometest = PowerPlant(working_fluid=fluid)
        eff = sometest.calculate_efficiency(210, 0.1, 65)
        if not np.isnan(eff):
            sometest.generate_diagram()
            sometest.plot_process(fn=fluid)
            sometest.print_result()
        else:
            print('+' * 75)
            print(fluid)
            print('+' * 75)
    except:
        print('+' * 75)
        print(fluid)
        print('+' * 75)
        pass
