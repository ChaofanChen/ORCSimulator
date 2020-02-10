# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 15:01:36 2019

@author: Chaofan
"""
import os, time, glob
import CoolProp
from CoolProp.Plots import PropertyPlot
from CoolProp.CoolProp import PropsSI

def compute(working_type, Q_w, Q_s, T_b_p, T_b_i, p_b, x_c, T_env, eta):
    type_wf = working_type
    T_wf_after_evaporator = T_b_p - 50
    # working fluid temperature out of the condenser decided by the condenser capcbility
    T_wf_after_condenser = T_env + 9
    # energy obtained from brine; '0' means liquid, '1' represents vapor
    H_b_production = PropsSI('H', 'T', T_b_p, 'Q', 0, 'Water') * Q_w + PropsSI('H', 'T', T_b_p, 'Q', 1, 'Water') * Q_s
    H_b_injection = (Q_w + Q_s) * PropsSI('H', 'T', T_b_i, 'Q', 0, 'Water')
    delt_P_b = (H_b_production - H_b_injection) * eta

    H_wf_before_preheater = PropsSI('H', 'T', T_wf_after_condenser, 'Q', 0, type_wf)
    H_wf_after_evaporator = PropsSI('H', 'T', T_wf_after_evaporator, 'Q', 1, type_wf)
    # mass flow rate of the working fluid
    Q_wf = delt_P_b / (H_wf_after_evaporator - H_wf_before_preheater)
    # pressure of the working fluid
    p_wf_before_pump = PropsSI('P', 'T', T_wf_after_condenser, 'Q', 1, type_wf)
    p_wf_after_pump = PropsSI('P', 'T', T_wf_after_evaporator, 'Q', 1, type_wf)

    # print(p_wf_before_pump, p_wf_after_pump, Q_wf)
    # isentropic expansion in the turbine, further calculating outlet temperature of the turbine
    # entropy of inlet working fluid of the turbine
    S_wf_before_turbine = PropsSI('S', 'H', H_wf_after_evaporator, 'P', p_wf_after_pump, type_wf)
    # turbine outlet temperature
    T_wf_after_turbine = PropsSI('T', 'S', S_wf_before_turbine, 'P', p_wf_before_pump, type_wf)
    # enthalpy of the outlet working fluid from the turbine
    H_wf_after_turbine = PropsSI('H', 'T', T_wf_after_turbine, 'P', p_wf_before_pump, type_wf)
    # Considering efficiency of the turbine (80-85%) and generator (95%), the output electricity (kW) can be calculated.
    P_generator = Q_wf * (H_wf_after_evaporator - H_wf_after_turbine) * 0.80 * 0.95 / 1000
    # energy from brine
    P_brine = H_b_production - H_b_injection
    # efficiency of the geothermal power plant (%)
    eta_electricity = P_generator * 1000 / P_brine
    # show key temperatures and efficiency of the eletricity generation
    # print(T_wf_after_evaporator, T_wf_after_turbine, P_generator, eta_electricity)

    ts_plot = PropertyPlot(type_wf, 'Ts', tp_limits='ORC')
    ts_plot.calc_isolines(CoolProp.iQ, num=2)
    ts_plot.calc_isolines(CoolProp.iP, iso_range=[p_wf_before_pump/1000, p_wf_after_pump/1000], num=2, rounding=True)
    ts_plot.draw()
    ts_plot.props[CoolProp.iP]['color'] = 'green'
    ts_plot.props[CoolProp.iP]['lw'] = '1'
    ts_plot.title(r'$T,s$ Graph for working fluid')
    ts_plot.xlabel(r'$s$ [kJ/kg K]')
    ts_plot.ylabel(r'$T$ [K]')
    ts_plot.grid()

    if not os.path.isdir('static'):
        os.mkdir('static')
    else:
        # Remove old plot files
        for filename in glob.glob(os.path.join('static', '*.png')):
            os.remove(filename)
    # Use time since Jan 1, 1970 in filename in order make
    # a unique filename that the browser has not chached
    plotfile = os.path.join('static', str(time.time()) + '.png')
    ts_plot.savefig(plotfile)
    return plotfile, P_generator, T_wf_after_evaporator, T_wf_after_turbine, p_wf_before_pump / 100000, p_wf_after_pump / 100000, Q_wf, eta_electricity

if __name__ == '__main__':
    compute(1, 0.1, 1, 20)
