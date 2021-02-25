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
fluids = ['R245fa', 'R600', 'R245CA', 'R123', 'Isopentane', 'n-Pentane', 'R113', 'R141B', 'R11']
for fluid in fluids:
    print('+' * 75)
    PP = geothermal_orc_design.PowerPlant(working_fluid=fluid)
    PP.documented = False

    print('Working fluid:', fluid)
    state = CP.AbstractState('HEOS', fluid)
    T_crit = state.trivial_keyed_output(CP.iT_critical) - 273.15
    print('Critical temperature: {} °C'.format(round(T_crit, 4)))

    sensitivity_analysis_without_ihe = pd.DataFrame(columns=['turbine inlet pressure', 'power_output', 'thermal_efficiency', 'net_power', 'net_efficiency', 'T_i', 'Q_IHE', 'Q_Brine_EV', 'geo_steam_share'])

    Q_range = -np.linspace(3e7, 1e3, 10)
    x_range = np.linspace(0.05, 0.20, 7)

    PP.nw.get_comp('internal heat exchanger').set_attr(pr1=0.98, pr2=0.98)
    for x in x_range:
        PP.run_simulation(brine_evap_Td=-5, IHE_sizing=.999, T_air_hot=15, geo_steam_share=x)

        sensitivity_analysis_without_ihe.loc[len(sensitivity_analysis_without_ihe)] = [
            PP.get_p_before_turbine(), -PP.get_power(), PP.get_efficiency(),
            -PP.get_net_power(), PP.get_net_efficiency(),
            PP.get_T_reinjection(),
            PP.get_internal_heat_exchanger_heat(),
            PP.get_brine_evaporator_heat(),
            x
        ]

    print(sensitivity_analysis_without_ihe)

    geothermal_orc_design.plot_sensitivity_analysis(
        sensitivity_analysis_without_ihe,
        fn='with_working_fluid_of_' + fluid,
        x='geo_steam_share',
        x_label='Geosteam mass fraction',
        y1='net_power', y2='T_i',
        y1_label='Net power output in MW', y2_label='Brine reinjection temperature in °C')
