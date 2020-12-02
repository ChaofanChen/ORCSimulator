#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 14:12:16 2020

@author: chencha
"""
import CoolProp as CP
import pandas as pd
import numpy as np
import geothermal_orc_design_fwi_chaofan 

# -----different working fluids at the same design logic (comparative study)-----------------------------
fluids = ['R600'] # 'R600', 'R245fa', 'R245CA', 'R11', 'Isopentane', 'n-Pentane', 'R123', 'R141B', 'R113'
sensitivity_analysis_Q_ihe = pd.DataFrame(columns=['power_output', 'thermal_efficiency', 'T_i'])
for fluid in fluids:
    print('+' * 75)
    print('Working fluid:', fluid)
    state = CP.AbstractState('HEOS', fluid)
    T_crit = state.trivial_keyed_output(CP.iT_critical) - 273.15
    print('Critical temperature: {} °C'.format(round(T_crit, 4)))

    Q_ihes = np.linspace(11, 13, 10)#8879392.8237
    p_turs = np.linspace(0, 0, 1)
    Td_bp_con =2.5437
    ihe_ttd_ls=np.linspace(0, 10, 10)
    for p_tur in p_turs:
        PowerPlantWithIHE = geothermal_orc_design_fwi_chaofan.PowerPlant(working_fluid=fluid)
        eff = PowerPlantWithIHE.calculate_efficiency_with_ihe(140, 200, 0.1, 70, p_tur, 1, 10, 10)
        sensitivity_analysis_Q_ihe.loc[p_tur, 'power_output'], \
        sensitivity_analysis_Q_ihe.loc[p_tur, 'thermal_efficiency'], \
        sensitivity_analysis_Q_ihe.loc[p_tur, 'T_i']=PowerPlantWithIHE.print_result()
#        PowerPlantWithIHE.plot_Ts(fn=fluid, Td_bp_cond=2)
print(sensitivity_analysis_Q_ihe)
geothermal_orc_design_fwi_chaofan.plot_sensitivity_analysis(sensitivity_analysis_Q_ihe,
                                                            fn='with_working_fluid_of_' + fluid, 
                                                            kw='Q_ihe with working fluid of '+fluid+' [W]')