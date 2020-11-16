#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 14:35:25 2020

@author: chencha
"""

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
import logging


logger.define_logging(screen_level=logging.ERROR)

def plot_sensitivity_analysis(sensitivity_analysis, fn='fluid', kw='Td_bp'):
    fig, ax = plt.subplots()
    ax.plot(sensitivity_analysis.index, sensitivity_analysis['power_output'], color='blue', marker="o")
    ax.set(xlabel=kw +' [°C]', ylabel='Power output [MW]')
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
    
class PowerPlantWithoutIHE():

    def __init__(self, working_fluid):

        self.working_fluid = working_fluid
        fluids = ['water', self.working_fluid, 'air']
        self.nw = network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # geo parameters

        geo_mass_flow = 200
        geo_steam_share = 0.1
        T_brine_in = 140
        T_reinjection = 70

        # ambient parameters

        T_amb = 5
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
        feed_working_fluid_pump = pump('feed working fluid pump')
        geo_merge = merge('brine merge')

        tur = turbine('turbine')

        ls_valve = valve('live steam valve')

        # ihe = heat_exchanger('internal heat exchanger')

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
        tur_cond = connection(tur, 'out1', air_cond, 'in1', label='tur_cond')
        self.nw.add_conns(ls_in, lsv_tur, tur_cond)

        # condenser to steam generator
        cond_fwp = connection(air_cond, 'out1', feed_working_fluid_pump, 'in1', label='cond_fwp')
        fwp_eco = connection(feed_working_fluid_pump, 'out1', eco, 'in2', label='fwp_eco')
        self.nw.add_conns(cond_fwp, fwp_eco)

        # steam generator
        eco_dr = connection(eco, 'out2', dr, 'in1', label='eco_dr')
        dr_esp = connection(dr, 'out1', evap_splitter, 'in1')
        esp_eb = connection(evap_splitter, 'out1', evap_brine, 'in2')
        esp_es = connection(evap_splitter, 'out2', evap_steam, 'in2')
        eb_em = connection(evap_brine, 'out2', evap_merge, 'in1')
        es_em = connection(evap_steam, 'out2', evap_merge, 'in2')
        em_dr = connection(evap_merge, 'out1', dr, 'in2')
        ls_out = connection(dr, 'out2', orc_cc, 'in1')
        self.nw.add_conns(eco_dr, dr_esp, esp_eb, esp_es, eb_em, es_em, em_dr, ls_out)

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
        cond_fwp.set_attr(fluid={self.working_fluid: 1, 'air': 0, 'water': 0})
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
        tur_cond.set_attr(h0=ws_stable_h0)
        tur_cond.set_attr(Td_bp=8, design=['Td_bp'], p0=PSI('P', 'T', T_amb + 273.15, 'Q', 1, self.working_fluid) / 1e5)
        fwp_eco.set_attr(h=ref(cond_fwp, 1, 1e3))

        # steam generator
        gs_es.set_attr(m=geo_mass_flow * geo_steam_share, T=T_brine_in, x=1, p0=5)
        gb_eb.set_attr(m=geo_mass_flow * (1 - geo_steam_share),
                       p=PSI('P', 'T', T_brine_in + 273.15, 'Q', 1, 'water') / 1e5, T=T_brine_in, state='l')

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
        air_cond.set_attr(pr1=1, pr2=1, ttd_u=10, design=['ttd_u'], offdesign=['kA_char'])
        feed_working_fluid_pump.set_attr(design=['eta_s'], offdesign=['eta_s_char'])

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
        # tur_cond.set_attr(h=None)
        fwp_eco.set_attr(h=None)
        eb_gm.set_attr(T=None)
        tur_cond.set_attr(Td_bp=None)
        eco_dr.set_attr(Td_bp=None)
        
        self.nw.connections['geosteam'].set_attr(m=20, T=140)
        self.nw.connections['geobrine'].set_attr(m=180, T=140)
        self.nw.connections['tur_cond'].set_attr(Td_bp=2)
        self.nw.connections['eco_dr'].set_attr(Td_bp=-2)

    def calculate_efficiency(self, geo_mass_flow, geo_steam_fraction, T_reinjection, Td_bp_cond, Td_bp_pre):

        self.nw.connections['geosteam'].set_attr(m=geo_mass_flow * geo_steam_fraction)
        self.nw.connections['geobrine'].set_attr(m=geo_mass_flow * (1 - geo_steam_fraction))
#        self.nw.connections['reinjection'].set_attr(T=T_reinjection)
        self.nw.connections['tur_cond'].set_attr(Td_bp=Td_bp_cond)
        self.nw.connections['eco_dr'].set_attr(Td_bp=Td_bp_pre)
        self.nw.solve('design')
        self.nw.print_results()

        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            return np.nan
        else:
            return self.nw.busses['power output'].P.val / self.nw.busses['thermal input'].P.val

    def calculate_efficiency_opt(self, x):

        self.nw.connections['tur_cond'].set_attr(Td_bp=x[0])
        self.nw.connections['eco_dr'].set_attr(Td_bp=x[1])

        self.nw.solve('design')
        # self.nw.print_results()

        for cp in self.nw.components.values():
            if isinstance(cp, heat_exchanger):
                if cp.Q.val > 0:
                    return np.nan

        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            return np.nan
        else:
            return self.nw.busses['power output'].P.val / self.nw.busses['thermal input'].P.val
        
    def print_result(self):
        for cp in self.nw.components.values():
            if isinstance(cp, heat_exchanger):
                if cp.Q.val > 0:
                    return 0, 0

        eta_th = self.nw.busses['power output'].P.val / self.nw.busses['thermal input'].P.val
        power = self.nw.busses['power output'].P.val
        print('Power output: {} MW'.format(round(power / 1e6, 4)))
        print('Thermal efficiency: {} %'.format(round(eta_th * 100, 4)))
        return power / 1e6, eta_th * 100, self.nw.connections['reinjection'].T.val

    def plot_Ts(self, fn='fluid', Td_bp_cond=1):

        T_before_turbine = self.nw.connections['lsv_tur'].T.val
        s_before_turbine = PSI('S', 'T', T_before_turbine + 273.15, 'Q', 1, self.working_fluid)

        T_after_turbine = self.nw.connections['tur_cond'].T.val
        s_after_turbine = PSI('S', 'T', T_after_turbine + 273.15, 'P', self.nw.connections['tur_cond'].p.val * 1e5, self.working_fluid)

        T_after_condenser = self.nw.connections['cond_fwp'].T.val
        s_after_condenser = PSI('S', 'T', T_after_condenser + 273.15, 'Q', 0, self.working_fluid)

        T_after_pump = self.nw.connections['fwp_eco'].T.val
        s_after_pump = PSI('S', 'T', T_after_pump + 273.15, 'P', self.nw.connections['fwp_eco'].p.val * 1e5, self.working_fluid)

        T_after_preheater = self.nw.connections['eco_dr'].T.val
        s_after_preheater = PSI('S', 'T', T_after_preheater + 273.15, 'P', self.nw.connections['eco_dr'].p.val * 1e5, self.working_fluid)

        state = CP.AbstractState('HEOS', self.working_fluid)
        T_crit = state.trivial_keyed_output(CP.iT_critical)
        df = pd.DataFrame(columns=['s_l', 's_g', 's_iso_P0', 's_iso_P1', 's_iso_P_top', 's_iso_P_bottom'])

        P0 = self.nw.connections['cond_fwp'].p.val * 1e5
        P1 = self.nw.connections['lsv_tur'].p.val * 1e5
        T_range = np.linspace(273.15, T_crit, 1000)
        for T in T_range:
            df.loc[T, 's_l'] = PSI('S', 'T', T, 'Q', 0, self.working_fluid)
            df.loc[T, 's_g'] = PSI('S', 'T', T, 'Q', 1, self.working_fluid)
            df.loc[T, 's_iso_P0'] = PSI('S', 'T', T, 'P', P0, self.working_fluid)
            df.loc[T, 's_iso_P1'] = PSI('S', 'T', T, 'P', P1, self.working_fluid)

        T_range_evaporator = np.linspace(self.nw.connections['eco_dr'].T.val + 273.15, self.nw.connections['lsv_tur'].T.val + 273.15 + 0.1, 100)
        for T in T_range_evaporator:
            df.loc[T, 's_iso_P_top'] = PSI('S', 'T', T, 'P', P1, self.working_fluid)

        T_steam_wf_low_P = PSI('T', 'P', P0, 'Q', 1, self.working_fluid)
        s_steam_wf_low_P = PSI('S', 'P', P0, 'Q', 1, self.working_fluid)
        T_range_condenser = np.linspace(T_steam_wf_low_P + Td_bp_cond, self.nw.connections['cond_fwp'].T.val + 273.15 - 0.1, 200)
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

        Temp = [T_before_turbine, T_after_turbine, T_steam_wf_low_P - 273.15, T_after_condenser,
                T_after_preheater] # , T_after_pump
        entropy = [s_before_turbine, s_after_turbine, s_steam_wf_low_P, s_after_condenser, s_after_preheater] # , s_after_pump
        n = ['1', '2', '', '3', '4'] 

        ax.scatter(entropy, Temp, color='red', s=0.5)
        for i, txt in enumerate(n):
            ax.annotate(txt, (entropy[i], Temp[i]), fontsize=12, textcoords="offset points", xytext=(0,6), horizontalalignment='right')  # , ha='center'
        for i in range(0, 1, 1):
            plt.plot(entropy[i:i + 2], Temp[i:i + 2], 'ro-', lw=2)
        for i in range(2, 5, 1):
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
#        plt.show()

#fluids = ['R245CA'] # 'R600', 'R245fa', 'R245CA', 'R11', 'Isopentane', 'n-Pentane', 'R123', 'R141B', 'R113'
#Td_bp_conds = np.linspace(20, 30.3, 30)
#for fluid in fluids:
##    try:
#        # for some testing
#        print('+' * 75)
#        sen_analy_Td_bp_cond_without_IHE = pd.DataFrame(columns=['power_output', 'thermal_efficiency', 'T_i'])
#        WithoutIHE = PowerPlantWithoutIHE(working_fluid=fluid)
#        print('Working fluid:', fluid)
#        state = CP.AbstractState('HEOS', fluid)
#        T_crit = state.trivial_keyed_output(CP.iT_critical) - 273.15
#        print('Critical temperature: {} °C'.format(round(T_crit, 4)))
#        for Td_bp_cond in Td_bp_conds:
#            eff_without_IHE = WithoutIHE.calculate_efficiency(200, 0.1, 70, Td_bp_cond, -2)
#            sen_analy_Td_bp_cond_without_IHE.loc[Td_bp_cond, 'power_output'],\
#            sen_analy_Td_bp_cond_without_IHE.loc[Td_bp_cond, 'thermal_efficiency'],\
#            sen_analy_Td_bp_cond_without_IHE.loc[Td_bp_cond, 'T_i']=WithoutIHE.print_result()
#            WithoutIHE.plot_Ts(fn=fluid, Td_bp_cond=Td_bp_cond)
#        plot_sensitivity_analysis(sen_analy_Td_bp_cond_without_IHE, fn=fluid, kw='Td_bp_condense_without_IHE')
#        
#        print(sen_analy_Td_bp_cond_without_IHE)
##    except:
##        print('+' * 75)
##        print(fluid)
##        print('+' * 75)
##        pass
        
        
import pygmo as pg

class optimization_problem():

    def fitness(self, x):
        f1 = 1 / self.model.calculate_efficiency_opt(x)
        ci1 = -x[0] + x[1]
        print(x)
        return [f1, ci1]

    def get_nobj(self):
        """Return number of objectives."""
        return 1

    # equality constraints
    def get_nec(self):
        return 0

    # inequality constraints
    def get_nic(self):
        return 1

    def get_bounds(self):
        """Return bounds of decision variables."""
        return ([2, -5], [30, 0])

optimize = optimization_problem()
optimize.model = PowerPlantWithoutIHE(working_fluid='Isopentane')
prob = pg.problem(optimize)
num_gen = 15

pop = pg.population(prob, size=10)
algo = pg.algorithm(pg.ihs(gen=num_gen))

result = {'champion': [], 'efficiency': [], 'generation': [],
          'extraction 1': [], 'extraction 2': []}

for gen in range(num_gen):
    result["generation"].append(gen)
    result["champion"].append(100/pop.champion_f[0])

    decision_var = pop.get_x()
    for pressure in decision_var:
        result['extraction 1'].append(pressure[0])
        result['extraction 2'].append(pressure[1])

    fitness = pop.get_f()
    for efficiency in fitness:
        result['efficiency'].append(100/efficiency[0])

    print()
    print('Evolution: {}'.format(gen))
    print('Thermal Efficiency: {} %'.format(round(100 / pop.champion_f[0], 4)))
    pop = algo.evolve(pop)

print()
print('Thermal Efficiency: {} %'.format(round(100 / pop.champion_f[0], 4)))
print('Td_bp before condenser: {} °C'.format(round(pop.champion_x[0], 4)))
print('Td_bp after preheater: {} °C'.format(round(pop.champion_x[1], 4)))

# scatter plot
cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(result['extraction 2'], result['extraction 1'], linewidth=0.25,
                 c=result['efficiency'], cmap=cm, alpha=0.5, edgecolors='black')
plt.scatter(pop.champion_x[1], pop.champion_x[0], marker='x', linewidth=1,
            c='red')
plt.annotate('Optimum', xy=(pop.champion_x[1], pop.champion_x[0]),
             xytext=(pop.champion_x[1]+2, pop.champion_x[0]+2),
             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',
                             color='red')
             )
plt.ylabel('Td_bp before condenser in °C')
plt.xlabel('Td_bp after preheater in °C')
plt.colorbar(sc, label='Cycle thermal efficiency (%)')
plt.savefig("Td_bp_optimization.png")
plt.show()