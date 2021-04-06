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


def golden_ratio_search(function, get_param, a, b, tol=1e-5, direction='min', param_to_opt=None, func_params={}):
    """Golden ratio search"""

    invphi = (5 ** 0.5 - 1) / 2  # 1 / phi
    invphi2 = (3 - 5 ** 0.5) / 2  # 1 / phi^2

    (a, b) = (min(a, b), max(a, b))
    h = b - a
    if h <= tol:
        return np.array([a, b])

    # Required steps to achieve tolerance
    n = int(np.ceil(np.log(tol / h) / np.log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h

    if direction == 'min':
        func_params[param_to_opt] = c
        function(**func_params)
        yc = get_param()
        func_params[param_to_opt] = d
        function(**func_params)
        yd = get_param()
    else:
        func_params[param_to_opt] = c
        function(**func_params)
        yc = 1 / get_param()
        func_params[param_to_opt] = d
        function(**func_params)
        yd = 1 / get_param()

    for k in range(n-1):
        if yc < yd:
            b = d
            d = c
            yd = yc
            h = invphi * h
            c = a + invphi2 * h
            func_params[param_to_opt] = c
            if direction == 'min':
                function(**func_params)
                yc = get_param()
            else:
                function(**func_params)
                yc = 1 / get_param()
        else:
            a = c
            c = d
            yc = yd
            h = invphi * h
            d = a + invphi * h
            func_params[param_to_opt] = d
            if direction == 'min':
                function(**func_params)
                yd = get_param()
            else:
                function(**func_params)
                yd = 1 / get_param()

    if yc < yd:
        return np.array([a, d])
    else:
        return np.array([c, b])


fluids = ['R245CA', 'Isopentane', 'n-Pentane', 'R123', 'R113', 'R141B', 'R11', 'R245fa', 'R600']
optima = pd.DataFrame(index=fluids)

for fluid in fluids:
    print('+' * 75)
    ORC = geothermal_orc_design.PowerPlant(working_fluid=fluid)
    ORC.documented = False
    print('Working fluid:', fluid)
    state = CP.AbstractState('HEOS', fluid)
    T_crit = state.trivial_keyed_output(CP.iT_critical) - 273.15
    print('Critical temperature: {} °C'.format(round(T_crit, 4)))
    ORC.nw.get_comp('internal heat exchanger').set_attr(pr1=1, pr2=1)

    p_opt = golden_ratio_search(
        ORC.run_simulation, ORC.get_power,
        a=-3e7, b=0, tol=1, param_to_opt='Q_brine_ev',
        func_params={'Q_ihe': 0, 'T_air_hot': 15})
    optima.loc[fluid, r'$p_1$ in bar'] = ORC.get_T_before_turbine()
    optima.loc[fluid, r'$p_2$ in bar'] = ORC.get_p_after_turbine()
    optima.loc[fluid, r'$T_2$ in °C'] = ORC.get_T_after_turbine()
    optima.loc[fluid, r'$T_\mathrm{reinj}$ in °C'] = ORC.get_T_reinjection()
    optima.loc[fluid, r'$\dot{Q}_\mathrm{BEv}$ in MW'] = -ORC.get_brine_evaporator_heat() / 1e6
    optima.loc[fluid, r'$P$ in MW'] = -ORC.get_power() / 1e6
    optima.loc[fluid, r'$\eta_\mathrm{th}$ in \%'] = ORC.get_efficiency() * 100
    optima.loc[fluid, r'$P_\mathrm{el}$ in MW'] = -ORC.get_net_power() / 1e6
    optima.loc[fluid, r'$\eta_\mathrm{net}$ in MW'] = ORC.get_net_efficiency() * 100

print(optima.to_latex(escape=False, na_rep='-', float_format='%.3f'))
