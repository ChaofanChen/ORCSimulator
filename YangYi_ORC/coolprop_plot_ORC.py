import CoolProp as CP
from CoolProp.CoolProp import PropsSI
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

state = CP.AbstractState('HEOS', 'Isopentane')
T_crit = state.trivial_keyed_output(CP.iT_critical)
df = pd.DataFrame(columns=['s_l', 's_g', 's_iso_P0', 's_iso_P1'])
P0=1000000
P1=200000
T_range = np.geomspace(273.15, T_crit, 1000)
for T in T_range:
    df.loc[T, 's_l'] = PropsSI('S', 'T', T, 'Q', 0, 'Isopentane')
    df.loc[T, 's_g'] = PropsSI('S', 'T', T, 'Q', 1, 'Isopentane')
    df.loc[T, 's_iso_P0'] = PropsSI('S', 'T', T, 'P', P0, 'Isopentane')
    df.loc[T, 's_iso_P1'] = PropsSI('S', 'T', T, 'P', P1, 'Isopentane')

print(df)
fig, ax = plt.subplots()
ax.plot(df['s_g']/1000, df.index - 273.15, color='black')
ax.plot(df['s_l']/1000, df.index - 273.15, color='black')
ax.plot(df['s_iso_P0']/1000, df.index - 273.15, color='green')
ax.plot(df['s_iso_P1']/1000, df.index - 273.15, color='green')
ax.set(xlabel='Specific entropy [kJ/kg K]', ylabel='Temperature [K]',
       title='T,s Graph for working fluid')
ax.grid()
plt.savefig('ts_plot_new.png')
plt.show()

#
# for back_end in ['HEOS', 'IF97']:
#
#     state = CP.AbstractState(back_end, 'water')
#
#     h = 146644.8016353955
#     p = 5628.620143029655
#     state.update(CP.HmassP_INPUTS, h, p)
#     print(state.smass())
#     # p = 5628.621143029655
#     # state.update(CP.HmassP_INPUTS, h, p)
#     # print(state.smass())
#
#     p_min = 5600
#     p_max = 5700
#
#     p_range = np.geomspace(p_min, p_max, 10)
#
#     df = pd.DataFrame(columns=['s', 'T', 'rho'])
#
#     for p in p_range:
#         try:
#             state.update(CP.HmassP_INPUTS, h, p)
#             df.loc[p, 's'] = state.smass()
#             df.loc[p, 'T'] = state.T()
#             df.loc[p, 'rho'] = state.rhomass()
#         except ValueError:
#             df.loc[p, 's'] = np.nan
#             df.loc[p, 'T'] = np.nan
#             df.loc[p, 'rho'] = np.nan
#
#     fig, ax = plt.subplots(3, 1)
#
#     ax[0].scatter(df.index, df['s'], marker='x')
#     ax[0].grid()
#     ax[0].set_xlabel('pressure in Pa')
#     ax[0].set_ylabel('entropy in J/kgK')
#     ax[1].scatter(df.index, df['T'], marker='x')
#     ax[1].grid()
#     ax[1].set_xlabel('pressure in Pa')
#     ax[1].set_ylabel('temperature in K')
#     ax[2].scatter(df.index, df['rho'], marker='x')
#     ax[2].grid()
#     ax[2].set_xlabel('pressure in Pa')
#     ax[2].set_ylabel('rho in kg/m3')
#     plt.tight_layout()
#     fig.savefig('diagram_' + back_end +'.pdf')
#     fig.savefig('diagram_' + back_end +'.png')



# import CoolProp
# from CoolProp.Plots import PropertyPlot
# from CoolProp.Plots import SimpleRankineCycle
# pp = PropertyPlot('Isopentane', 'Ts', unit_system='EUR')
# pp.calc_isolines(CoolProp.iQ, num=2)
# cycle = SimpleRankineCycle('Isopentane', 'TS', unit_system='EUR')
# T0 = 27.7+273.15
# pp.state.update(CoolProp.QT_INPUTS,0.0,T0-15)
# p0 = pp.state.keyed_output(CoolProp.iP)
# T2 = 117.3+273.15
# pp.state.update(CoolProp.QT_INPUTS,1.0,T2+10)
# p2 = pp.state.keyed_output(CoolProp.iP)
# pp.calc_isolines(CoolProp.iP, [p0/1e5,p2/1e5], num=2, rounding=True)
# pp.props[CoolProp.iP]['color'] = 'green'
# pp.props[CoolProp.iP]['lw'] = '1'
# cycle.simple_solve(T0, p0, T2, p2, 1, 0.4, SI=True)
# cycle.steps = 50
# sc = cycle.get_state_changes()
# pp.draw_process(sc, line_opts={'color':'red', 'lw':1.5})
# import matplotlib.pyplot as plt
# plt.close(cycle.figure)
# pp.show()





