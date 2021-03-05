from tespy.connections import Connection, Bus, Ref
from tespy.networks import Network
from tespy.components import (
    HeatExchanger, Pump, Turbine, Source, Sink, CycleCloser, Splitter,
    Merge, Condenser, Drum, Valve, HeatExchangerSimple, Compressor
)
from tespy.tools.helpers import UserDefinedEquation
from tespy.tools.fluid_properties import h_mix_pT, T_mix_ph
from tespy.tools import logger

#from fluprodia import FluidPropertyDiagram
import CoolProp as CP
from CoolProp.CoolProp import PropsSI as PSI
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import math
import logging
import tespy


logger.define_logging(screen_level=logging.WARNING)


def plot_sensitivity_analysis(
        sensitivity_analysis, fn='fluid', y1='', y2='',
        y1_label='', y2_label='', x_label=''):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    ax.plot(sensitivity_analysis.index, sensitivity_analysis[y1], color='blue', marker="o")
    ax.set(xlabel=x_label, ylabel=y1_label)
    ax2=ax.twinx()
    ax2.plot(sensitivity_analysis.index, sensitivity_analysis[y2], color='red', marker="*")
    ax2.set(ylabel=y2_label)
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
    plt.savefig('diff_' + x_label + '_plot_' + fn + '.png')


def desuperheat(ude):

    return (
        ude.conns[0].h.val_SI - ude.conns[1].h.val_SI -
        ude.params['distance'] * (
            ude.conns[0].h.val_SI - h_mix_pT(
                ude.conns[1].get_flow(),
                T_mix_ph(ude.conns[2].get_flow()))))


def desuperheat_deriv(ude):
    ude.jacobian[ude.conns[0]][2] = 1 - ude.params['distance']
    ude.jacobian[ude.conns[1]][1] = ude.numeric_deriv('p', 1)
    ude.jacobian[ude.conns[1]][2] = -1
    ude.jacobian[ude.conns[2]][1] = ude.numeric_deriv('p', 2)
    ude.jacobian[ude.conns[2]][2] = ude.numeric_deriv('h', 2)
    return ude.jacobian


class PowerPlant():

    def __init__(self, working_fluid):

        self.working_fluid = working_fluid
        fluids = ['water', self.working_fluid, 'air']
        self.nw = Network(fluids=fluids)
        self.nw.set_attr(p_unit='bar', T_unit='C', h_unit='kJ / kg')

        # geo parameters

        geo_mass_flow = 200
        geo_steam_share = 0.1
        T_brine_in = 140
#        T_reinjection = 70

        # ambient parameters

        self.T_amb = 5
        self.p_amb = 0.6

        # main components

        geo_steam = Source('steam source')
        geo_brine = Source('brine source')
        geo_reinjection = Sink('reinjection well')

        air_in = Source('ambient air source')
        air_out = Sink('ambient air sink')
        air_fan = Compressor('ambient air fan')
        air_cond = Condenser('main condenser')

        orc_cc = CycleCloser('orc cycle closer')

        # evap_steam = Condenser('steam evaporator')
        evap_splitter = Splitter('splitter evaporation')
        evap_merge = Merge('merge evaporation')
        evap_steam = Condenser('steam evaporator')
        # geo_brine_pump = Pump('geobrine pump')
        evap_brine = HeatExchanger('brine evaporator')
        dr = Drum('drum')

        eco = HeatExchanger('economiser')
        feed_working_fluid_pump = Pump('feed working fluid pump')
        geo_merge = Merge('brine merge')

        tur = Turbine('turbine')

        ihe = HeatExchanger('internal heat exchanger')

        # busses
        net_power = Bus('net power output')
        net_power.add_comps(
            {'comp': tur, 'char': 0.97},
            {'comp': feed_working_fluid_pump, 'char': 0.97, 'base': 'bus'},
            # {'comp': geo_brine_pump, 'char': 0.97, 'base': 'bus'},
            {'comp': air_fan, 'char': 0.97, 'base': 'bus'}
        )

        ORC_power_bus = Bus('ORC power')
        ORC_power_bus.add_comps(
            {'comp': tur}, {'comp': feed_working_fluid_pump}
        )

        geothermal_bus = Bus('thermal input')
        geothermal_bus.add_comps(
            {'comp': eco, 'char': -1}, {'comp': evap_brine, 'char': -1},
            {'comp': evap_steam, 'char': -1}
        )

        self.nw.add_busses(net_power, ORC_power_bus, geothermal_bus)

        # turbine to condenser
        ls_in = Connection(orc_cc, 'out1', tur, 'in1', label='lsv_tur')
        tur_ihe = Connection(tur, 'out1', ihe, 'in1', label='tur_ihe')
        ihe_cond = Connection(ihe, 'out1', air_cond, 'in1', label='ihe_cond')
        self.nw.add_conns(ls_in, tur_ihe, ihe_cond)

        # condenser to steam generator
        cond_fwp = Connection(air_cond, 'out1', feed_working_fluid_pump, 'in1', label='cond_fwp')
        fwp_ihe = Connection(feed_working_fluid_pump, 'out1', ihe, 'in2', label='fwp_ihe')
        self.nw.add_conns(cond_fwp, fwp_ihe)

        # steam generator
        ihe_eco = Connection(ihe, 'out2', eco, 'in2', label='ihe_eco')
        eco_dr = Connection(eco, 'out2', dr, 'in1', label='eco_dr')
        dr_esp = Connection(dr, 'out1', evap_splitter, 'in1')
        esp_eb = Connection(evap_splitter, 'out1', evap_brine, 'in2')
        esp_es = Connection(evap_splitter, 'out2', evap_steam, 'in2')
        eb_em = Connection(evap_brine, 'out2', evap_merge, 'in1', label='eb_em')
        es_em = Connection(evap_steam, 'out2', evap_merge, 'in2', label='es_em')
        em_dr = Connection(evap_merge, 'out1', dr, 'in2')
        ls_out = Connection(dr, 'out2', orc_cc, 'in1')
        self.nw.add_conns(ihe_eco, eco_dr, dr_esp, esp_eb, esp_es, eb_em, es_em, em_dr, ls_out)

        # air cold side
        air_in_fan = Connection(air_in, 'out1', air_fan, 'in1')
        fan_cond = Connection(air_fan, 'out1', air_cond, 'in2')
        cond_air_hot = Connection(air_cond, 'out2', air_out, 'in1', 'air hot')
        self.nw.add_conns(air_in_fan, fan_cond, cond_air_hot)

        # geo source
        gs_es = Connection(geo_steam, 'out1', evap_steam, 'in1', label='geosteam')
        es_gm = Connection(evap_steam, 'out1',  geo_merge, 'in1')
        gb_gm = Connection(geo_brine, 'out1', geo_merge, 'in2', label='geobrine')
        gm_eb = Connection(geo_merge, 'out1', evap_brine, 'in1', label='geobrine mix')
        self.nw.add_conns(gs_es, es_gm, gb_gm, gm_eb)

        eb_eco = Connection(evap_brine, 'out1', eco, 'in1', 'brine to eco')
        eco_gr = Connection(eco, 'out1', geo_reinjection, 'in1', label='reinjection')
        self.nw.add_conns(eb_eco, eco_gr)

        # fluid settings
        ihe_eco.set_attr(fluid={self.working_fluid: 1.0, 'air': 0.0, 'water': 0.0})
        air_in_fan.set_attr(fluid={self.working_fluid: 0.0, 'air': 1.0, 'water': 0.0})
        gs_es.set_attr(fluid={self.working_fluid: 0.0, 'air': 0.0, 'water': 1.0})
        gb_gm.set_attr(fluid={self.working_fluid: 0.0, 'air': 0.0, 'water': 1.0})

        # connection parameters
        ls_stable_p0 = PSI('P', 'T', T_brine_in + 273.15, 'Q', 1, self.working_fluid) / 1e5
        ls_in.set_attr(p0=ls_stable_p0)
        ws_stable_h0 = (
            PSI('H', 'T', self.T_amb + 273.15, 'Q', 1, self.working_fluid) + 0.5 * (
                PSI('H', 'T', T_brine_in + 273.15, 'Q', 1, self.working_fluid) -
                PSI('H', 'T', self.T_amb + 273.15, 'Q', 1, self.working_fluid)
            )
        ) / 1e3
        tur_ihe.set_attr(h=ws_stable_h0)
        ihe_cond.set_attr(Td_bp=8, design=['Td_bp'], p0=PSI('P', 'T', self.T_amb + 273.15, 'Q', 1, self.working_fluid) / 1e5)
        fwp_ihe.set_attr(h=Ref(cond_fwp, 1, 1))

        # steam generator
        gs_es.set_attr(m=geo_mass_flow * geo_steam_share, T=T_brine_in, x=1, p0=5)
        gb_gm.set_attr(m=geo_mass_flow * (1 - geo_steam_share), T=T_brine_in, x=0)

        em_dr.set_attr()
        eb_em.set_attr(x=0.5)
        es_em.set_attr(x=0.5, design=['x'])
        eb_eco.set_attr(h=Ref(gm_eb, 1, -50))

        eco_dr.set_attr(Td_bp=-2)
        # eco_dr.set_attr(x=0)

        # main condenser
        air_in_fan.set_attr(p=self.p_amb, T=self.T_amb)
        cond_air_hot.set_attr(T=self.T_amb + 15, p=self.p_amb)

        # component parameters
        # turbines
        tur.set_attr(design=['eta_s'], offdesign=['cone', 'eta_s_char'])
        # condensing
        ihe.set_attr(pr1=0.98, pr2=0.98, offdesign=['kA_char'])
        air_cond.set_attr(pr1=1, pr2=0.995, ttd_u=10, design=['ttd_u'], offdesign=['kA_char'])
        air_fan.set_attr(eta_s=0.6)
        feed_working_fluid_pump.set_attr(design=['eta_s']) # has already set eta_s, thus drop 'offdesign=['eta_s_char']'

        # steam generator
        evap_steam.set_attr(offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        evap_brine.set_attr(pr1=0.98, ttd_l=8, offdesign=['kA_char'])  # no pr2 due to drum pressure balance
        eco.set_attr(pr1=0.98, pr2=0.98)
        # geo_brine_pump.set_attr(eta_s=0.75, design=['eta_s'], offdesign=['eta_s_char'])

        self.nw.set_attr(iterinfo=False)
        self.nw.solve('design')
        self.nw.save('stable_' + self.working_fluid)
        # self.nw.print_results()
        tur.set_attr(eta_s=0.9)
        feed_working_fluid_pump.set_attr(eta_s=0.75)
        tur_ihe.set_attr(h=None)
        fwp_ihe.set_attr(h=None)
        eb_eco.set_attr(h=None, T=Ref(gm_eb, 1, -10))
        self.nw.get_conn('ihe_cond').set_attr(Td_bp=None)

        self.ude_IHE_size = UserDefinedEquation(
            label='ihe deshuperheat ratio',
            func=desuperheat, deriv=desuperheat_deriv,
            conns=[
                self.nw.get_conn('tur_ihe'),
                self.nw.get_conn('ihe_cond'),
                self.nw.get_conn('fwp_ihe')],
            params={'distance': 0.0}
        )

    def run_simulation(
        self, p_before_tur=None, Q_ihe=None, Q_brine_ev=None,
        T_reinjection=None, brine_evap_Td=None, T_air_hot=None, IHE_sizing=None, 
        geo_mass=200, geo_steam_share=0.1):
        self.nw.get_comp('internal heat exchanger').set_attr(Q=Q_ihe)
        self.nw.get_conn('lsv_tur').set_attr(p=p_before_tur)
        self.nw.get_conn('reinjection').set_attr(T=T_reinjection)
        self.nw.get_comp('brine evaporator').set_attr(Q=Q_brine_ev)
        self.nw.get_conn('geosteam').set_attr(m=geo_mass*geo_steam_share)
        self.nw.get_conn('geobrine').set_attr(m=geo_mass*(1-geo_steam_share))

        if brine_evap_Td is not None:
            self.nw.get_conn('brine to eco').set_attr(T=Ref(self.nw.get_conn('geobrine mix'), 1, brine_evap_Td))
        else:
            self.nw.get_conn('brine to eco').set_attr(T=None)

        if T_air_hot is not None:
            self.nw.get_conn('air hot').set_attr(T=self.T_amb + T_air_hot)
        else:
            self.nw.get_conn('air hot').set_attr(T=None)

        if IHE_sizing is None:
            if self.ude_IHE_size in self.nw.user_defined_eq.values():
                self.nw.del_ude(self.ude_IHE_size)
        else:
            if self.ude_IHE_size not in self.nw.user_defined_eq.values():
                self.nw.add_ude(self.ude_IHE_size)
            self.ude_IHE_size.params['distance'] = IHE_sizing

        try:
            self.nw.solve('design')
#            self.nw.print_results()
        except ValueError:
            self.nw.res = [1]
            pass

    def check_simulation(self, value):
        if self.nw.lin_dep or self.nw.res[-1] > 1e-3:
            self.nw.solve('design', init_path='stable_' + self.working_fluid, init_only=True)
            return np.nan
        else:
            for cp in self.nw.comps['object']:
                if isinstance(cp, HeatExchanger):
                    if cp.Q.val > 0:
                        return np.nan
                    elif cp.kA.val <= 0:
                        return np.nan
                    elif cp.ttd_l.val <= 0:
                        return np.nan
        return value

    def get_power(self):
        return self.check_simulation(self.nw.busses['ORC power'].P.val)

    def get_net_power(self):
        return self.check_simulation(self.nw.busses['net power output'].P.val)

    def get_efficiency(self):
        return self.check_simulation(-self.nw.busses['ORC power'].P.val / self.nw.busses['thermal input'].P.val)

    def get_net_efficiency(self):
        return self.check_simulation(-self.nw.busses['net power output'].P.val / self.nw.busses['thermal input'].P.val)

    def get_T_reinjection(self):
        return self.check_simulation(self.nw.get_conn('reinjection').T.val)

    def get_p_before_turbine(self):
        return self.check_simulation(self.nw.get_conn('lsv_tur').p.val)

    def get_p_after_turbine(self):
        return self.check_simulation(self.nw.get_conn('tur_ihe').p.val)

    def get_T_after_turbine(self):
        return self.check_simulation(self.nw.get_conn('tur_ihe').T.val)

    def get_brine_evaporator_heat(self):
        return self.check_simulation(self.nw.get_comp('brine evaporator').Q.val)

    def get_internal_heat_exchanger_heat(self):
        return self.check_simulation(self.nw.get_comp('internal heat exchanger').Q.val)

    def get_working_fluid_mass_flow(self):
        return self.check_simulation(self.nw.get_conn('tur_ihe').m.val)

    def plot_Ts(self, fn='fluid', Td_bp_cond=1):

        T_before_turbine = self.nw.get_conn('lsv_tur').T.val
        s_before_turbine = PSI('S', 'T', T_before_turbine + 273.15, 'Q', 1, self.working_fluid)

        T_after_turbine = self.nw.get_conn('tur_ihe').T.val
        s_after_turbine = PSI('S', 'T', T_after_turbine + 273.15, 'P', self.nw.get_conn('tur_ihe').p.val * 1e5, self.working_fluid)

        T_before_condenser = self.nw.get_conn('ihe_cond').T.val
        s_before_condenser = PSI('S', 'T', T_before_condenser + 273.15, 'P', self.nw.get_conn('ihe_cond').p.val * 1e5, self.working_fluid)

        T_after_condenser = self.nw.get_conn('cond_fwp').T.val
        s_after_condenser = PSI('S', 'T', T_after_condenser + 273.15, 'Q', 0, self.working_fluid)

        T_after_pump = self.nw.get_conn('fwp_ihe').T.val
        s_after_pump = PSI('S', 'T', T_after_pump + 273.15, 'P', self.nw.get_conn('fwp_ihe').p.val * 1e5, self.working_fluid)

        T_after_ihe = self.nw.get_conn('ihe_eco').T.val
        s_after_ihe = PSI('S', 'T', T_after_ihe + 273.15, 'P', self.nw.get_conn('ihe_eco').p.val * 1e5, self.working_fluid)

        T_after_preheater = self.nw.get_conn('eco_dr').T.val
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

        # s_brine_in = PSI('S', 'T', self.nw.get_conn('geobrine').T.val + 273.15, 'P',
        #                  self.nw.get_conn('geobrine').p.val, 'water')
        # s_brine_out = PSI('S', 'T', self.nw.get_conn('reinjection').T.val + 273.15, 'P',
        #                   self.nw.get_conn('reinjection').p.val * 100000, 'water')
        # ax.scatter(s_brine_in, self.nw.get_conn('geobrine').T.val)
        # ax.scatter(s_brine_out, self.nw.get_conn('reinjection').T.val)

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
