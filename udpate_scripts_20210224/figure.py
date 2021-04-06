#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:25:36 2020

@author: chencha
"""

import pandas as pd
import matplotlib.pyplot as plt

fluids = ['R245fa', 'R600', 'R245CA', 'R123', 'Isopentane', 'n-Pentane', 'R113', 'R141B', 'R11'] #, 'Isobutane',

for fluid in fluids:
    df = pd.read_csv('diff_T_before_tur_' + fluid + '_without_ihe.csv', 
                         header=0, names=['T_before_tur', 'power_output','thermal_efficiency',
                                          'net_power', 'net_efficiency', 'T_i', 'Q_IHE', 'Q_Brine_EV'])
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    ax.plot(df['T_before_tur'], df['power_output']/1e6, color='blue', marker="o") #, label = 'Gross power output'
#    ax.plot(df['p_before_tur'], df['net_power']/1e6, color='blue', marker="x")
    ax.set(xlabel= 'Turbine inlet temperature with ' + fluid + ' (°C)', ylabel='Gross power output amount (MW)')
#    plt.ylim(14, 19)
    ax2=ax.twinx()
    ax2.plot(df['T_before_tur'], df['T_i'], color='red', marker="*", label = 'Re-injection temperature')
    ax2.set(ylabel='Re-injection temperature (°C)')
    plt.ylim(45, 100)
    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(18)
    ax.xaxis.label.set_size(18)
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=15)
    # ax.set_ticklabel(exclude_overlapping=True)
    ax.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_color('red')
    ax2.yaxis.label.set_size(18)
    ax2.tick_params(axis="y", labelsize=15)
    ax2.tick_params(axis='y', colors='red')
    ax.grid()
#    ax.legend()
#    ax2.legend()
#    plt.show()
    # fig.autofmt_xdate()
    plt.savefig('diff_T_before_tur_' + fluid + '.pdf')
#    print(df)