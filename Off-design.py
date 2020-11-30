#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 13:55:59 2020

@author: chencha
"""
import pandas as pd
import numpy as np
import geothermal_orc_design_fwi_chaofan 

# ------------off-design performance analysis----------------------------------------------------------
off_design_performance_f = pd.DataFrame(columns=['power_output', 'thermal_efficiency', 'T_i'])
off_design_performance_T_pro = pd.DataFrame(columns=['power_output', 'thermal_efficiency', 'T_i'])
steam_mass_fractions = np.linspace(0.1, 0.09, 10)
T_productions = np.linspace(142, 138, 11)
PowerPlantOffDesign = geothermal_orc_design_fwi_chaofan.PowerPlant(working_fluid='Isopentane')
 
# for steam_mass_fraction in steam_mass_fractions:
#     offeff = PowerPlantOffDesign.calculate_efficiency_off_design(140, 200, steam_mass_fraction, 70)
#     off_design_performance_f.loc[steam_mass_fraction, 'power_output'],\
#     off_design_performance_f.loc[steam_mass_fraction, 'thermal_efficiency']= PowerPlantOffDesign.print_result()
# plot_sensitivity_analysis(off_design_performance_f, fn='Isopentane', kw='Steam mass fraction')

for T_production in T_productions:
    offeff = PowerPlantOffDesign.calculate_efficiency_off_design(T_production, 200, 0.1, 70)
    off_design_performance_T_pro.loc[T_production, 'power_output'], \
    off_design_performance_T_pro.loc[T_production, 'thermal_efficiency'], \
    off_design_performance_T_pro.loc[T_production, 'T_i']= PowerPlantOffDesign.print_result()

geothermal_orc_design_fwi_chaofan.plot_sensitivity_analysis(off_design_performance_T_pro, fn='Isopentane', kw='T_pro')