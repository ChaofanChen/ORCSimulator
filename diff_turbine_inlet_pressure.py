#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 14:05:02 2020

@author: chencha
"""
import CoolProp as CP
import pandas as pd
import numpy as np
import geothermal_orc_design_fwi_chaofan 

# --------sensitivity analysis for every single parameter--------------------------
fluids = ['R600'] #, 'Isobutane','R245fa', 'R600', 'R245CA', 'R123', 'Isopentane', 'n-Pentane', 'R113', 'R141B', 'R11'
for fluid in fluids:
    print('+' * 75)
    PowerPlantWithIHE = geothermal_orc_design_fwi_chaofan.PowerPlant(working_fluid=fluid)
    print('Working fluid:', fluid)
    state = CP.AbstractState('HEOS', fluid)
    T_crit = state.trivial_keyed_output(CP.iT_critical) - 273.15
    print('Critical temperature: {} °C'.format(round(T_crit, 4)))

    sensitivity_analysis_without_ihe = pd.DataFrame(columns=['power_output', 'thermal_efficiency', 'T_i'])
    p_before_turs = np.linspace(24, 26.3, 10)
    for p_before_tur in p_before_turs:
        eff = PowerPlantWithIHE.calculate_efficiency_without_ihe(p_before_tur)
        sensitivity_analysis_without_ihe.loc[p_before_tur, 'power_output'], \
        sensitivity_analysis_without_ihe.loc[p_before_tur, 'thermal_efficiency'],\
        sensitivity_analysis_without_ihe.loc[p_before_tur, 'T_i']=PowerPlantWithIHE.print_result()
    print(sensitivity_analysis_without_ihe)
#    sensitivity_analysis_without_ihe.to_csv('diff_p_before_tur_' + fluid + '_without_ihe.csv')
    geothermal_orc_design_fwi_chaofan.plot_sensitivity_analysis(sensitivity_analysis_without_ihe, 
                                                                fn=fluid, kw='without_ihe')
