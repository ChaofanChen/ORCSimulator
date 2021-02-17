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


invphi = (5 ** 0.5 - 1) / 2  # 1 / phi
invphi2 = (3 - 5 ** 0.5) / 2  # 1 / phi^2

def golden_ratio_search(function, a, b, tol=1e-5, direction='min'):
    """Golden ratio search"""

    (a, b) = (min(a, b), max(a, b))
    h = b - a
    if h <= tol:
        return np.array([a, b])

    # Required steps to achieve tolerance
    n = int(np.ceil(np.log(tol / h) / np.log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h

    if direction == 'min':
        yc = function(c)
        yd = function(d)
    else:
        yc = 1 / function(c)
        yd = 1 / function(d)

    for k in range(n-1):
        if yc < yd:
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            if direction == 'min':
                yc = function(c)
            else:
                yc = 1 / function(c)
        else:
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            if direction == 'min':
                yd = function(d)
            else:
                yd = 1 / function(d)

    if yc < yd:
        return np.array([a, d])
    else:
        return np.array([c, b])

# --------sensitivity analysis for every single parameter--------------------------
fluids = ['R245fa', 'R600', 'R245CA', 'R123', 'Isopentane', 'n-Pentane', 'R113', 'R141B', 'R11']
optima = pd.DataFrame(index=fluids)

for fluid in fluids:
    print('+' * 75)
    PowerPlantWithIHE = geothermal_orc_design_fwi_chaofan.PowerPlant(working_fluid=fluid)
    PowerPlantWithIHE.documented = False
    print('Working fluid:', fluid)
    state = CP.AbstractState('HEOS', fluid)
    T_crit = state.trivial_keyed_output(CP.iT_critical) - 273.15
    print('Critical temperature: {} °C'.format(round(T_crit, 4)))

    sensitivity_analysis_without_ihe = pd.DataFrame(columns=['power_output', 'thermal_efficiency', 'T_i'])
    state.update(
        CP.QT_INPUTS, 1,
        PowerPlantWithIHE.nw.get_conn('geobrine').T.val_SI -
        30 - PowerPlantWithIHE.nw.get_comp('brine evaporator').ttd_l.val)
    p_min = state.p() / 1e5
    state.update(
        CP.QT_INPUTS, 1,
        PowerPlantWithIHE.nw.get_conn('geobrine').T.val_SI -
        0.1 - PowerPlantWithIHE.nw.get_comp('brine evaporator').ttd_l.val)
    p_max = state.p() / 1e5
    p_before_turs = np.linspace(p_min, p_max, 10)

    p_opt = golden_ratio_search(PowerPlantWithIHE.calculate_power_without_ihe, p_min, p_max, 1e-5, direction='max')
    optima.loc[fluid, r'$p_1$ in bar'] = PowerPlantWithIHE.get_p_before_turbine()
    optima.loc[fluid, r'$p_2$ in bar'] = PowerPlantWithIHE.get_p_after_turbine()
    optima.loc[fluid, r'$T_2$ in °C'] = PowerPlantWithIHE.get_T_after_turbine()
    optima.loc[fluid, r'$T_\mathrm{reinj}$ in °C'] = PowerPlantWithIHE.get_T_reinjection()
    optima.loc[fluid, r'$\dot{Q}_\mathrm{BEv}$ in MW'] = PowerPlantWithIHE.get_brine_evaporator_heat() / 1e6
    optima.loc[fluid, r'$\eta_\mathrm{th}$ in \%'] = PowerPlantWithIHE.get_efficiency() * 100
    optima.loc[fluid, r'$P$ in MW'] = PowerPlantWithIHE.get_power() / 1e6

print(optima.to_latex(escape=False, na_rep='-', float_format='%.3f'))
