from tespy.connections import connection, bus, ref
from tespy.tools import char_line
from tespy.networks import network
from tespy.components import (
    heat_exchanger, pump, turbine, source, sink, cycle_closer, splitter,
    merge, condenser, drum, valve, heat_exchanger_simple
)
from tespy.tools import logger

#from fluprodia import FluidPropertyDiagram
import CoolProp as CP
from CoolProp.CoolProp import PropsSI as PSI
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import math
import logging


logger.define_logging(screen_level=logging.ERROR)

def plot_sensitivity_analysis(sensitivity_analysis, fn='fluid', kw='Td_bp'):
    fig, ax = plt.subplots()
    ax.plot(sensitivity_analysis.index, sensitivity_analysis['power_output'], color='blue', marker="o")
    ax.set(xlabel=kw, ylabel='Net power output [MW]')
    ax2=ax.twinx()
    ax2.plot(sensitivity_analysis.index, sensitivity_analysis['thermal_efficiency'], color='red', marker="*")
    ax2.set(ylabel='Thermal efficiency [%]')
    # plt.ylim(10, 20)
    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(12)
    ax.xaxis.label.set_size(12)
    # ax.set_ticklabel(exclude_overlapping=True)
    ax.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_color('red')
    ax2.yaxis.label.set_size(12)
    ax2.tick_params(axis='y', colors='red')
    ax.grid()
    # plt.show()
    # fig.autofmt_xdate()
    plt.savefig('diff_' + kw + '_plot_' + fn + '.png')

class PowerPlant():

    def __init__(self, working_fluid):

        self.working_fluid = working_fluid
        fluids = ['water', self.working_fluid, 'air']
        self.nw = network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # geo parameters

        geo_mass_flow = 200
        geo_steam_share = 0.1
        T_brine_in = 140
#        T_reinjection = 70

        # ambient parameters

        T_amb = 5
        p_amb = 0.6

        # main components

        geo_steam = source('steam source')
        geo_brine = source('brine source')
        geo_reinjection = sink('reinjection well')

        air_in = source('ambient air source')
        air_out = sink('ambient air sink')
        air_cond = condenser('main condenser')

        orc_cc = cycle_closer('orc cycle closer')

        # evap_steam = condenser('steam evaporator')
        evap_splitter = splitter('splitter evaporation')
        evap_merge = merge('merge evaporation')
        evap_steam = condenser('steam evaporator')
        geo_steam_pump = pump('geosteam condensate pump')
        evap_brine = heat_exchanger('brine evaporator')
        dr = drum('drum')

        eco = heat_exchanger('economiser')
        feed_working_fluid_pump = pump('feed working fluid pump')
        geo_merge = merge('brine merge')

        tur = turbine('turbine')

        ls_valve = valve('live steam valve')

        ihe = heat_exchanger('internal heat exchanger')

        # busses
        power_bus = bus('power output')
        power_bus.add_comps(
            {'c': tur, 'char': -1},
            {'c': feed_working_fluid_pump, 'char': -1}, {'c': geo_steam_pump, 'char': -1}
        )

        geothermal_bus = bus('thermal input')
        geothermal_bus.add_comps(
            {'c': eco, 'char': -1}, {'c': evap_brine, 'char': -1},
            {'c': evap_steam, 'char': -1}
        )

        self.nw.add_busses(power_bus, geothermal_bus)

        # turbine to condenser
        ls_in = connection(orc_cc, 'out1', ls_valve, 'in1')
        lsv_tur = connection(ls_valve, 'out1', tur, 'in1', label='lsv_tur')
        tur_ihe = connection(tur, 'out1', ihe, 'in1', label='tur_ihe')
        ihe_cond = connection(ihe, 'out1', air_cond, 'in1', label='ihe_cond')
        self.nw.add_conns(ls_in, lsv_tur, tur_ihe, ihe_cond)

        # condenser to steam generator
        cond_fwp = connection(air_cond, 'out1', feed_working_fluid_pump, 'in1', label='cond_fwp')
        fwp_ihe = connection(feed_working_fluid_pump, 'out1', ihe, 'in2', label='fwp_ihe')
        self.nw.add_conns(cond_fwp, fwp_ihe)

        # steam generator
        ihe_eco = connection(ihe, 'out2', eco, 'in2', label='ihe_eco')
        eco_dr = connection(eco, 'out2', dr, 'in1', label='eco_dr')
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
        ihe_cond.set_attr(Td_bp=8, design=['Td_bp'], p0=PSI('P', 'T', T_amb + 273.15, 'Q', 1, self.working_fluid) / 1e5)
        fwp_ihe.set_attr(h=ref(cond_fwp, 1, 1e3))

        # steam generator
        gs_es.set_attr(m=geo_mass_flow * geo_steam_share, T=T_brine_in, x=1, p0=5)
        gb_eb.set_attr(m=geo_mass_flow * (1 - geo_steam_share), T=T_brine_in, x=0)

        em_dr.set_attr()
        eb_em.set_attr(x=0.4)
        es_em.set_attr(x=0.6, design=['x'])
        eb_gm.set_attr(T=T_brine_in - 20)

#        eco_dr.set_attr(Td_bp=-2)
        eco_dr.set_attr(x=0)

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
        feed_working_fluid_pump.set_attr(design=['eta_s']) # has already set eta_s, thus drop 'offdesign=['eta_s_char']'

        # steam generator
        evap_steam.set_attr(pr1=0.99, offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        evap_brine.set_attr(pr1=1, ttd_l=10, offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        eco.set_attr(pr1=1, pr2=1)
        geo_steam_pump.set_attr(eta_s=0.75, design=['eta_s'], offdesign=['eta_s_char'])

        self.nw.set_attr(iterinfo=False)
        self.nw.solve('design')
#        self.nw.print_results()
        tur.set_attr(eta_s=0.9)
        feed_working_fluid_pump.set_attr(eta_s=0.75)
        tur_ihe.set_attr(h=None)
        fwp_ihe.set_attr(h=None)
        eb_gm.set_attr(T=None)

    
    def calculate_efficiency_opt_without_ihe(self, x, working_fluid):

        self.working_fluid = working_fluid
        fluids = ['water', self.working_fluid, 'air']
        self.nw = network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # geo parameters

        geo_mass_flow = 200
        geo_steam_share = 0.1
        T_brine_in = 140
#        T_reinjection = 70

        # ambient parameters

        T_amb = 5
        p_amb = 0.6

        # main components

        geo_steam = source('steam source')
        geo_brine = source('brine source')
        geo_reinjection = sink('reinjection well')

        air_in = source('ambient air source')
        air_out = sink('ambient air sink')
        air_cond = condenser('main condenser')

        orc_cc = cycle_closer('orc cycle closer')

        # evap_steam = condenser('steam evaporator')
        evap_splitter = splitter('splitter evaporation')
        evap_merge = merge('merge evaporation')
        evap_steam = condenser('steam evaporator')
        geo_steam_pump = pump('geosteam condensate pump')
        evap_brine = heat_exchanger('brine evaporator')
        dr = drum('drum')

        eco = heat_exchanger('economiser')
        feed_working_fluid_pump = pump('feed working fluid pump')
        geo_merge = merge('brine merge')

        tur = turbine('turbine')

        ls_valve = valve('live steam valve')

        ihe = heat_exchanger('internal heat exchanger')

        # busses
        power_bus = bus('power output')
        power_bus.add_comps(
            {'c': tur, 'char': -1},
            {'c': feed_working_fluid_pump, 'char': -1}, {'c': geo_steam_pump, 'char': -1}
        )

        geothermal_bus = bus('thermal input')
        geothermal_bus.add_comps(
            {'c': eco, 'char': -1}, {'c': evap_brine, 'char': -1},
            {'c': evap_steam, 'char': -1}
        )

        self.nw.add_busses(power_bus, geothermal_bus)

        # turbine to condenser
        ls_in = connection(orc_cc, 'out1', ls_valve, 'in1')
        lsv_tur = connection(ls_valve, 'out1', tur, 'in1', label='lsv_tur')
        tur_ihe = connection(tur, 'out1', ihe, 'in1', label='tur_ihe')
        ihe_cond = connection(ihe, 'out1', air_cond, 'in1', label='ihe_cond')
        self.nw.add_conns(ls_in, lsv_tur, tur_ihe, ihe_cond)

        # condenser to steam generator
        cond_fwp = connection(air_cond, 'out1', feed_working_fluid_pump, 'in1', label='cond_fwp')
        fwp_ihe = connection(feed_working_fluid_pump, 'out1', ihe, 'in2', label='fwp_ihe')
        self.nw.add_conns(cond_fwp, fwp_ihe)

        # steam generator
        ihe_eco = connection(ihe, 'out2', eco, 'in2', label='ihe_eco')
        eco_dr = connection(eco, 'out2', dr, 'in1', label='eco_dr')
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
        ihe_cond.set_attr(Td_bp=8, design=['Td_bp'], p0=PSI('P', 'T', T_amb + 273.15, 'Q', 1, self.working_fluid) / 1e5)
        fwp_ihe.set_attr(h=ref(cond_fwp, 1, 1e3))

        # steam generator
        gs_es.set_attr(m=geo_mass_flow * geo_steam_share, T=T_brine_in, x=1, p0=5)
        gb_eb.set_attr(m=geo_mass_flow * (1 - geo_steam_share), T=T_brine_in, x=0)

        em_dr.set_attr()
        eb_em.set_attr(x=0.5)
        es_em.set_attr(x=0.5, design=['x'])
        eb_gm.set_attr(T=T_brine_in - 20)

#        eco_dr.set_attr(Td_bp=-2)
        eco_dr.set_attr(x=0)

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
        feed_working_fluid_pump.set_attr(design=['eta_s']) # has already set eta_s, thus drop 'offdesign=['eta_s_char']'

        # steam generator
        evap_steam.set_attr(pr1=0.99, offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        evap_brine.set_attr(pr1=1, ttd_l=10, offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        eco.set_attr(pr1=1, pr2=1)
        geo_steam_pump.set_attr(eta_s=0.75, design=['eta_s'], offdesign=['eta_s_char'])

        self.nw.set_attr(iterinfo=False)
        self.nw.solve('design')
        # self.nw.print_results()
        tur.set_attr(eta_s=0.9)
        feed_working_fluid_pump.set_attr(eta_s=0.75)
        tur_ihe.set_attr(h=None)
        fwp_ihe.set_attr(h=None)
        eb_gm.set_attr(T=None)
        
        self.nw.components['internal heat exchanger'].set_attr(Q=0)
#        print(x[0], x[1])
        self.nw.connections['ihe_cond'].set_attr(Td_bp=None)
#        self.nw.connections['eco_dr'].set_attr(Td_bp=-0.01)
        self.nw.connections['lsv_tur'].set_attr(p=x[0])

        self.nw.solve('design')
#        self.nw.print_results()

#        for cp in self.nw.components.values():
#            if isinstance(cp, heat_exchanger):
#                if cp.Q.val > 0:
#                    return np.nan

        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            return np.nan
        else:
            return self.nw.busses['power output'].P.val / 1e6 #/ self.nw.busses['thermal input'].P.val

    def calculate_efficiency_off_design(self, T_production, geo_mass_flow, geo_steam_fraction, T_reinjection):

        self.nw.connections['geosteam'].set_attr(m=20, T=140)
        self.nw.connections['geobrine'].set_attr(m=180, T=140)
        self.nw.connections['reinjection'].set_attr(T=70)
        self.nw.connections['ihe_cond'].set_attr(Td_bp=2)
#        self.nw.connections['eco_dr'].set_attr(Td_bp=-2)
        self.nw.components['main condenser'].set_attr(ttd_u=10)
        self.nw.components['brine evaporator'].set_attr(ttd_l=10)
        self.nw.solve('design')
        self.nw.save('ORC')

        self.nw.connections['geosteam'].set_attr(m=geo_mass_flow * geo_steam_fraction, T=T_production)
        self.nw.connections['geobrine'].set_attr(m=geo_mass_flow * (1 - geo_steam_fraction), T=T_production)
        self.nw.connections['reinjection'].set_attr(T=T_reinjection)
        self.nw.solve('offdesign', design_path='ORC')
        # self.nw.print_results()

        # for cp in self.nw.components.values():
        #     if isinstance(cp, heat_exchanger):
        #         if cp.Q.val > 0:
        #             return np.nan
        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            return np.nan
        else:
            return self.nw.busses['power output'].P.val / self.nw.busses['thermal input'].P.val

    def calculate_efficiency_without_ihe(self, p_before_tur):

        self.nw.connections['ihe_cond'].set_attr(Td_bp=None)
#        self.nw.connections['eco_dr'].set_attr(Td_bp=-0.01)
        self.nw.components['internal heat exchanger'].set_attr(ttd_l=5)#Q=-1.1e7
        self.nw.connections['lsv_tur'].set_attr(p=p_before_tur)

        self.nw.solve('design')
        self.nw.print_results()

        for cp in self.nw.components.values():
            if isinstance(cp, heat_exchanger):
                if cp.Q.val > 0:
                    return np.nan

        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            return np.nan
        else:
            return self.nw.busses['power output'].P.val / self.nw.busses['thermal input'].P.val

    def calculate_efficiency_opt_with_ihe(self, x, working_fluid):

        self.working_fluid = working_fluid
        fluids = ['water', self.working_fluid, 'air']
        self.nw = network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # geo parameters

        geo_mass_flow = 200
        geo_steam_share = 0.1
        T_brine_in = 140
#        T_reinjection = 70

        # ambient parameters

        T_amb = 5
        p_amb = 0.6

        # main components

        geo_steam = source('steam source')
        geo_brine = source('brine source')
        geo_reinjection = sink('reinjection well')

        air_in = source('ambient air source')
        air_out = sink('ambient air sink')
        air_cond = condenser('main condenser')

        orc_cc = cycle_closer('orc cycle closer')

        # evap_steam = condenser('steam evaporator')
        evap_splitter = splitter('splitter evaporation')
        evap_merge = merge('merge evaporation')
        evap_steam = condenser('steam evaporator')
        geo_steam_pump = pump('geosteam condensate pump')
        evap_brine = heat_exchanger('brine evaporator')
        dr = drum('drum')

        eco = heat_exchanger('economiser')
        feed_working_fluid_pump = pump('feed working fluid pump')
        geo_merge = merge('brine merge')

        tur = turbine('turbine')

        ls_valve = valve('live steam valve')

        ihe = heat_exchanger('internal heat exchanger')

        # busses
        power_bus = bus('power output')
        power_bus.add_comps(
            {'c': tur, 'char': -1},
            {'c': feed_working_fluid_pump, 'char': -1}, {'c': geo_steam_pump, 'char': -1}
        )

        geothermal_bus = bus('thermal input')
        geothermal_bus.add_comps(
            {'c': eco, 'char': -1}, {'c': evap_brine, 'char': -1},
            {'c': evap_steam, 'char': -1}
        )

        self.nw.add_busses(power_bus, geothermal_bus)

        # turbine to condenser
        ls_in = connection(orc_cc, 'out1', ls_valve, 'in1')
        lsv_tur = connection(ls_valve, 'out1', tur, 'in1', label='lsv_tur')
        tur_ihe = connection(tur, 'out1', ihe, 'in1', label='tur_ihe')
        ihe_cond = connection(ihe, 'out1', air_cond, 'in1', label='ihe_cond')
        self.nw.add_conns(ls_in, lsv_tur, tur_ihe, ihe_cond)

        # condenser to steam generator
        cond_fwp = connection(air_cond, 'out1', feed_working_fluid_pump, 'in1', label='cond_fwp')
        fwp_ihe = connection(feed_working_fluid_pump, 'out1', ihe, 'in2', label='fwp_ihe')
        self.nw.add_conns(cond_fwp, fwp_ihe)

        # steam generator
        ihe_eco = connection(ihe, 'out2', eco, 'in2', label='ihe_eco')
        eco_dr = connection(eco, 'out2', dr, 'in1', label='eco_dr')
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
        ihe_cond.set_attr(Td_bp=8, design=['Td_bp'], p0=PSI('P', 'T', T_amb + 273.15, 'Q', 1, self.working_fluid) / 1e5)
        fwp_ihe.set_attr(h=ref(cond_fwp, 1, 1e3))

        # steam generator
        gs_es.set_attr(m=geo_mass_flow * geo_steam_share, T=T_brine_in, x=1, p0=5)
        gb_eb.set_attr(m=geo_mass_flow * (1 - geo_steam_share), T=T_brine_in, x=0)

        em_dr.set_attr()
        eb_em.set_attr(x=0.5)
        es_em.set_attr(x=0.5, design=['x'])
        eb_gm.set_attr(T=T_brine_in - 20)

#        eco_dr.set_attr(Td_bp=-2)
        eco_dr.set_attr(x=0)

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
        feed_working_fluid_pump.set_attr(design=['eta_s']) # has already set eta_s, thus drop 'offdesign=['eta_s_char']'

        # steam generator
        evap_steam.set_attr(pr1=0.99, offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        evap_brine.set_attr(pr1=1, ttd_l=10, offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        eco.set_attr(pr1=1, pr2=1)
        geo_steam_pump.set_attr(eta_s=0.75, design=['eta_s'], offdesign=['eta_s_char'])

        self.nw.set_attr(iterinfo=False)
        self.nw.solve('design')
        # self.nw.print_results()
        tur.set_attr(eta_s=0.9)
        feed_working_fluid_pump.set_attr(eta_s=0.75)
        tur_ihe.set_attr(h=None)
        fwp_ihe.set_attr(h=None)
        eb_gm.set_attr(T=None)
        
        self.nw.connections['lsv_tur'].set_attr(p=x[0])
        self.nw.components['internal heat exchanger'].set_attr(Q=x[1]*1e6)
        self.nw.connections['ihe_cond'].set_attr(Td_bp=None)
#        self.nw.connections['eco_dr'].set_attr(Td_bp=-0.01)

        self.nw.solve('design')
#        self.nw.print_results()

        for cp in self.nw.components.values():
            if isinstance(cp, heat_exchanger):
                if cp.Q.val > 0 or math.isnan(cp.ttd_l.val) < 0:
                    return np.nan

        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            return np.nan
        else:
            return self.nw.busses['power output'].P.val / 1e6 #/ self.nw.busses['thermal input'].P.val

    def calculate_efficiency_with_ihe(self, T_production, geo_mass_flow, geo_steam_fraction, T_reinjection, p_tur, Q_ihe, ttd_u_cond, ttd_l_evap):

        self.nw.components['internal heat exchanger'].set_attr(Q=0)
#        self.nw.connections['geosteam'].set_attr(m=geo_mass_flow * geo_steam_fraction, T=T_production)
#        self.nw.connections['geobrine'].set_attr(m=geo_mass_flow * (1 - geo_steam_fraction), T=T_production)
#        self.nw.connections['reinjection'].set_attr(T=57)
        self.nw.connections['ihe_cond'].set_attr(Td_bp=None)
#        self.nw.connections['eco_dr'].set_attr(Td_bp=Td_bp_eco)
#        self.nw.components['main condenser'].set_attr(ttd_u=ttd_u_cond)
#        self.nw.components['brine evaporator'].set_attr(ttd_l=ttd_l_evap)
        self.nw.connections['lsv_tur'].set_attr(p=26.3)

        self.nw.solve('design')
        self.nw.print_results()

        for cp in self.nw.components.values():
            if isinstance(cp, heat_exchanger):
                if cp.Q.val > 0 or math.isnan(cp.ttd_l.val) < 0:
                    return np.nan

        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            return np.nan
        else:
            return self.nw.busses['power output'].P.val / self.nw.busses['thermal input'].P.val

    
    def print_result(self):

        for cp in self.nw.components.values():
            if isinstance(cp, heat_exchanger):
                if cp.Q.val > 0:                 # or math.isnan(cp.kA.val)
                    return 0, 0, 0

        eta_th = self.nw.busses['power output'].P.val / self.nw.busses['thermal input'].P.val
        power = self.nw.busses['power output'].P.val
        print('Power output: {} MW'.format(round(power / 1e6, 4)))
        print('Thermal efficiency: {} %'.format(round(eta_th * 100, 4)))
        return power / 1e6, eta_th * 100, self.nw.connections['reinjection'].T.val
    
    def plot_Ts(self, fn='fluid', Td_bp_cond=1):

        T_before_turbine = self.nw.connections['lsv_tur'].T.val
        s_before_turbine = PSI('S', 'T', T_before_turbine + 273.15, 'Q', 1, self.working_fluid)

        T_after_turbine = self.nw.connections['tur_ihe'].T.val
        s_after_turbine = PSI('S', 'T', T_after_turbine + 273.15, 'P', self.nw.connections['tur_ihe'].p.val * 1e5, self.working_fluid)

        T_before_condenser = self.nw.connections['ihe_cond'].T.val
        s_before_condenser = PSI('S', 'T', T_before_condenser + 273.15, 'P', self.nw.connections['ihe_cond'].p.val * 1e5, self.working_fluid)

        T_after_condenser = self.nw.connections['cond_fwp'].T.val
        s_after_condenser = PSI('S', 'T', T_after_condenser + 273.15, 'Q', 0, self.working_fluid)

        T_after_pump = self.nw.connections['fwp_ihe'].T.val
        s_after_pump = PSI('S', 'T', T_after_pump + 273.15, 'P', self.nw.connections['fwp_ihe'].p.val * 1e5, self.working_fluid)

        T_after_ihe = self.nw.connections['ihe_eco'].T.val
        s_after_ihe = PSI('S', 'T', T_after_ihe + 273.15, 'P', self.nw.connections['ihe_eco'].p.val * 1e5, self.working_fluid)

        T_after_preheater = self.nw.connections['eco_dr'].T.val
        s_after_preheater = PSI('S', 'T', T_after_preheater + 273.15, 'Q', 0, self.working_fluid)

        state = CP.AbstractState('HEOS', self.working_fluid)
        T_crit = state.trivial_keyed_output(CP.iT_critical)
        df = pd.DataFrame(columns=['s_l', 's_g', 's_iso_P0', 's_iso_P1', 's_iso_P_top', 's_iso_P_bottom'])

        P0 = PSI('P', 'T', T_after_condenser + 273.15, 'Q', 0, self.working_fluid)
        P1 = PSI('P', 'T', T_before_turbine + 273.15, 'Q', 1, self.working_fluid)
        T_range = np.linspace(273.15, T_crit, 1000)
        for T in T_range:
            df.loc[T, 's_l'] = PSI('S', 'T', T, 'Q', 0, self.working_fluid)
            df.loc[T, 's_g'] = PSI('S', 'T', T, 'Q', 1, self.working_fluid)
            df.loc[T, 's_iso_P0'] = PSI('S', 'T', T, 'P', P0, self.working_fluid)
            df.loc[T, 's_iso_P1'] = PSI('S', 'T', T, 'P', P1, self.working_fluid)

        T_range_evaporator = np.linspace(T_after_preheater + 273.15 - 0.5, T_before_turbine + 273.15 + 0.5, 100)
        for T in T_range_evaporator:
            df.loc[T, 's_iso_P_top'] = PSI('S', 'T', T, 'P', P1, self.working_fluid)

        T_steam_wf_low_P = PSI('T', 'P', P0, 'Q', 1, self.working_fluid)
        s_steam_wf_low_P = PSI('S', 'P', P0, 'Q', 1, self.working_fluid)
        T_range_condenser = np.linspace(T_steam_wf_low_P + Td_bp_cond, T_after_condenser + 273.15 - 0.1, 200)
        for T in T_range_condenser:
            df.loc[T, 's_iso_P_bottom'] = PSI('S', 'T', T, 'P', P0, self.working_fluid)
        # print(df)

        fig, ax = plt.subplots()
        ax.plot(df['s_g'], df.index - 273.15, color='black')
        ax.plot(df['s_l'], df.index - 273.15, color='black')
        ax.plot(df['s_iso_P0'], df.index - 273.15, color='green')
        ax.plot(df['s_iso_P1'], df.index - 273.15, color='green')
        ax.plot(df['s_iso_P_top'], df.index - 273.15, color='red')
        ax.plot(df['s_iso_P_bottom'], df.index - 273.15, color='red')

        Temp = [T_before_turbine, T_after_turbine, T_before_condenser, T_steam_wf_low_P - 273.15, T_after_condenser,
                T_after_ihe, T_after_preheater] # , T_after_pump
        entropy = [s_before_turbine, s_after_turbine, s_before_condenser, s_steam_wf_low_P, s_after_condenser, s_after_ihe,
                   s_after_preheater] # , s_after_pump
        n = ['1', '2', '3', ' ', '4', '5', '6'] # , '7'

        ax.scatter(entropy, Temp, color='red', s=0.5)
        for i, txt in enumerate(n):
            ax.annotate(txt, (entropy[i], Temp[i]), fontsize=12, textcoords="offset points", xytext=(0,6), horizontalalignment='right')  # , ha='center'
        for i in range(0, 2, 1):
            plt.plot(entropy[i:i + 2], Temp[i:i + 2], 'ro-', lw=2)
        for i in range(4, 6, 1):
            plt.plot(entropy[i:i + 2], Temp[i:i + 2], 'ro-', lw=2)

        # s_brine_in = PSI('S', 'T', self.nw.connections['geobrine'].T.val + 273.15, 'P',
        #                  self.nw.connections['geobrine'].p.val, 'water')
        # s_brine_out = PSI('S', 'T', self.nw.connections['reinjection'].T.val + 273.15, 'P',
        #                   self.nw.connections['reinjection'].p.val * 100000, 'water')
        # ax.scatter(s_brine_in, self.nw.connections['geobrine'].T.val)
        # ax.scatter(s_brine_out, self.nw.connections['reinjection'].T.val)

        ax.set(xlabel='Specific entropy [J/kg K]', ylabel='Temperature [°C]')
        ax.grid()
        plt.savefig('ORC_Ts_plot_' + fn + '.png')
        # plt.show()

        
    # def plot_process(self, fn='somename'):
    #
    #     result_dict = {
    #         prop: [
    #             conn.get_attr(prop).val for conn in self.nw.conns.index
    #             if conn.fluid.val[self.working_fluid] == 1
    #         ] for prop in ['p', 'h', 's', 'T']
    #     }
    #
    #     self.diagram.set_limits(
    #         x_min=min(result_dict['h']) - 50,
    #         x_max=max(result_dict['h']) + 50,
    #         y_min=min(result_dict['p']) / 2,
    #         y_max=max(result_dict['p']) * 10
    #     )
    #     self.diagram.draw_isolines('logph')
    #     self.diagram.ax.scatter(result_dict['h'], result_dict['p'], zorder=100)
    #     self.diagram.save(fn + '_logph.pdf')
    #
    #     self.diagram.set_limits(
    #         x_min=min(result_dict['s']) - 50,
    #         x_max=max(result_dict['s']) + 50,
    #         y_min=min(result_dict['T']) - 25,
    #         y_max=PSI('T_critical', self.working_fluid) + 25 - 273.15
    #     )
    #     self.diagram.draw_isolines('Ts')
    #     self.diagram.ax.scatter(result_dict['s'], result_dict['T'], zorder=100)
    #     self.diagram.save(fn + '_Ts.pdf')
    #
    # def generate_diagram(self):
    #
    #     self.diagram = FluidPropertyDiagram(self.working_fluid)
    #     self.diagram.set_unit_system(T='°C', p='bar', h='kJ/kg')
    #     iso_T = np.arange(0, 201, 25)
    #     self.diagram.set_isolines(T=iso_T)
    #     self.diagram.calc_isolines()
    
    