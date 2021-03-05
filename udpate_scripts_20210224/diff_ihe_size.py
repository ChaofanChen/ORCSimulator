#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 14:05:02 2020

@author: chencha
"""
import CoolProp as CP
import pandas as pd
import numpy as np
import geothermal_orc_design
from tespy.tools.helpers import UserDefinedEquation
from tespy.tools.fluid_properties import h_mix_pT, T_mix_ph
import matplotlib.pyplot as plt

def desuperheat(ude):

    return (
        ude.conns[0].h.val_SI - ude.conns[1].h.val_SI -
        ude.params['distance'] * (
            ude.conns[0].h.val_SI - h_mix_pT(ude.conns[1].get_flow(), T_mix_ph(ude.conns[2].get_flow()))))


def desuperheat_deriv(ude):
    ude.jacobian[ude.conns[0]][2] = 1 - ude.params['distance']
    ude.jacobian[ude.conns[1]][1] = ude.numeric_deriv('p', 1)
    ude.jacobian[ude.conns[1]][2] = -1
    ude.jacobian[ude.conns[2]][1] = ude.numeric_deriv('p', 2)
    ude.jacobian[ude.conns[2]][2] = ude.numeric_deriv('h', 2)
    return ude.jacobian

# --------sensitivity analysis for every single parameter--------------------------
fluids = ['R245fa', 'R600'] #'R245fa', 'R600', 'R245CA', 'R123', 'Isopentane', 'n-Pentane', 'R113', 'R141B', 'R11'
for fluid in fluids:
    print('+' * 75)
    PP = geothermal_orc_design.PowerPlant(working_fluid=fluid)
    PP.documented = False

    # ihe_desuperheat_ude = UserDefinedEquation(
    #     label='ihe deshuperheat ratio',
    #     func=desuperheat, deriv=desuperheat_deriv,
    #     conns=[
    #         PP.nw.get_conn('tur_ihe'),
    #         PP.nw.get_conn('ihe_cond'),
    #         PP.nw.get_conn('fwp_ihe')],
    # # specify to 0.99999 for maximum physically possible heat extraction
    #     params={'distance': 0.0}
    # )
    #
    # PP.nw.add_ude(ihe_desuperheat_ude)
    print('Working fluid:', fluid)
    state = CP.AbstractState('HEOS', fluid)
    T_crit = state.trivial_keyed_output(CP.iT_critical) - 273.15
    print('Critical temperature: {} °C'.format(round(T_crit, 4)))

    sensitivity_analysis_without_ihe = pd.DataFrame(columns=['power_output', 'thermal_efficiency', 'net_power', 'net_efficiency', 'T_i', 'Q_IHE', 'Q_Brine_EV'])

    Q_range = -np.linspace(8e6, 0, 20)

#    PP.nw.get_comp('internal heat exchanger').set_attr(pr1=1, pr2=1)
    for Q in Q_range:
        PP.run_simulation(Q_brine_ev=-0.611e6, Q_ihe=Q, T_air_hot=15)

        sensitivity_analysis_without_ihe.loc[-PP.get_internal_heat_exchanger_heat()] = [
            -PP.get_power(), PP.get_efficiency()*100, -PP.get_net_power(),
            PP.get_net_efficiency()*100, PP.get_T_reinjection(),
            PP.get_internal_heat_exchanger_heat(),
            PP.get_brine_evaporator_heat()
        ]

    print(sensitivity_analysis_without_ihe)
    
    fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
    ax.plot(sensitivity_analysis_without_ihe.index/1e6, sensitivity_analysis_without_ihe['thermal_efficiency'], color='blue', marker="o", label='Gross thermal efficiency')
    ax.plot(sensitivity_analysis_without_ihe.index/1e6, sensitivity_analysis_without_ihe['net_efficiency'], color='g', marker="x", label='Net thermal efficiency')
    ax.set(xlabel= 'Heat exchange rate of IHE with ' + fluid + ' (MW)', ylabel='Thermal efficiency (%)')
    plt.ylim(10, 18)
    ax2=ax.twinx()
    ax2.plot(sensitivity_analysis_without_ihe.index/1e6, sensitivity_analysis_without_ihe['T_i'], color='black', marker="*", label='Re-injection temperature')
    ax2.set(ylabel='Re-injection temperature (°C)')
    plt.ylim(60, 75)
    plt.xlim(0, 8)
#    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(15)
    ax.xaxis.label.set_size(15)
    # ax.set_ticklabel(exclude_overlapping=True)
#    ax.tick_params(axis='y', colors='blue')
#    ax2.yaxis.label.set_color('red')
    ax2.yaxis.label.set_size(15)
#    ax2.tick_params(axis='y', colors='red')
    ax2.plot(np.linspace(0, 8, 5), (70,70,70,70,70), dashes=[2, 2], linewidth=3, color='red')
    ax.grid()
    ax.legend(loc='upper left')
    ax2.legend(loc='lower right')
#    plt.show()
    # fig.autofmt_xdate()
    plt.savefig('diff_IHE_size_with_' + fluid + '.png')
    
#    sensitivity_analysis_without_ihe.to_csv('diff_ihe_' + fluid + '.csv')

#    geothermal_orc_design.plot_sensitivity_analysis(
#        sensitivity_analysis_without_ihe,
#        fn='with_working_fluid_of_' + fluid,
#        y1='thermal_efficiency', y2='T_i',
#        y1_label='Net thermal efficiency (%)', y2_label='Re-injection temperature (°C)', x_label='Heat exchange rate of IHE (MW)')
