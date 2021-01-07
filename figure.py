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
    df = pd.read_csv('diff_p_before_tur_' + fluid + '_without_ihe.csv', 
                         header=0, names=['p_before_tur', 'power_output','thermal_efficiency', 'T_i'])
    fig, ax = plt.subplots()
    ax.plot(df['p_before_tur'], df['power_output'], color='blue', marker="o")
    ax.set(xlabel= 'Turbine inlet pressure with ' + fluid + ' [bar]', ylabel='Net power output [MW]')
#    plt.ylim(14, 19)
    ax2=ax.twinx()
    ax2.plot(df['p_before_tur'], df['T_i'], color='red', marker="*")
    ax2.set(ylabel='Re-injection temperature [°C]')
    plt.ylim(50, 95)
    ax.yaxis.label.set_color('blue')
    ax.yaxis.label.set_size(15)
    ax.xaxis.label.set_size(15)
    # ax.set_ticklabel(exclude_overlapping=True)
    ax.tick_params(axis='y', colors='blue')
    ax2.yaxis.label.set_color('red')
    ax2.yaxis.label.set_size(15)
    ax2.tick_params(axis='y', colors='red')
    ax.grid()
#    plt.show()
    # fig.autofmt_xdate()
    plt.savefig('diff_p_before_tur_' + fluid + '.pdf')
#    print(df)