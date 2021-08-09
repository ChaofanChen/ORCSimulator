#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 14:25:36 2020

@author: chencha
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

fluids = ['R245fa', 'R600', 'R245ca', 'Isopentane']

R245fa = pd.read_csv('IHE_sizing_with_R245fa.csv', 
                     header=0, names=['T_1', 'T_2','p_2',
                                      'T_35', 'Q_internal heat exchanger', 'gross power output',
                                      'net power output', 'thermal efficiency', 'net efficiency', 'IHE sizing factor'])
R600 = pd.read_csv('IHE_sizing_with_R600.csv', 
                     header=0, names=['T_1', 'T_2','p_2',
                                      'T_35', 'Q_internal heat exchanger', 'gross power output',
                                      'net power output', 'thermal efficiency', 'net efficiency', 'IHE sizing factor'])
R245ca = pd.read_csv('IHE_sizing_with_R245ca.csv', 
                     header=0, names=['T_1', 'T_2','p_2',
                                      'T_35', 'Q_internal heat exchanger', 'gross power output',
                                      'net power output', 'thermal efficiency', 'net efficiency', 'IHE sizing factor'])
Isopentane = pd.read_csv('IHE_sizing_with_Isopentane.csv', 
                     header=0, names=['T_1', 'T_2','p_2',
                                      'T_35', 'Q_internal heat exchanger', 'gross power output',
                                      'net power output', 'thermal efficiency', 'net efficiency', 'IHE sizing factor'])

fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
ax.plot(R245fa['Q_internal heat exchanger'], R245fa['net power output'], color='red', marker="o", label='R245fa')
ax.plot(R600['Q_internal heat exchanger'], R600['net power output'], color='blue', marker="x", label='R600')
ax.plot(R245ca['Q_internal heat exchanger'], R245ca['net power output'], color='black', marker="*", label='R245ca')
ax.plot(Isopentane['Q_internal heat exchanger'], Isopentane['net power output'], color='navy', marker="v", label='R601a (Isopentane)')
ax.set(xlabel= 'IHE heat exchange rate (MW)', ylabel='Net power output (MW)')
#plt.ylim(8, 14)
#ax2=ax.twinx()
#ax2.plot(R245fa['Q_internal heat exchanger'], R245fa['T_35'], color='black', marker="d", label='Re-injection temperature')
#ax2.set(ylabel='Re-injection temperature (°C)')
plt.ylim(10.5, 13)
plt.xlim(0, 12)
#ax.yaxis.label.set_color('blue')
ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18)
ax.tick_params(axis="x", labelsize=15)
ax.tick_params(axis="y", labelsize=15)
#ax.tick_params(axis='y', colors='blue')
#ax2.yaxis.label.set_size(18)
#ax2.tick_params(axis="y", labelsize=15)
#ax2.plot(np.linspace(0, 8, 5), (70,70,70,70,70), dashes=[2, 2], linewidth=3, color='red')
ax.grid()
ax.legend(frameon=False, prop={'size': 14})
plt.savefig('IHE_influence_net_power_output.pdf')


fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
ax.plot(R245fa['Q_internal heat exchanger'], R245fa['T_35'], color='red', marker="o", label='R245fa')
ax.plot(R600['Q_internal heat exchanger'], R600['T_35'], color='blue', marker="x", label='R600')
ax.plot(R245ca['Q_internal heat exchanger'], R245ca['T_35'], color='black', marker="*", label='R245ca')
ax.plot(Isopentane['Q_internal heat exchanger'], Isopentane['T_35'], color='navy', marker="v", label='R601a (Isopentane)')
ax.set(xlabel= 'IHE heat exchange rate (MW)', ylabel='Re-injection temperature (°C)')
ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18)
ax.tick_params(axis="x", labelsize=15)
ax.tick_params(axis="y", labelsize=15)
ax.yaxis.label.set_size(18)
plt.xlim(0, 12)
plt.ylim(60, 85)
ax.tick_params(axis="y", labelsize=15)
ax.plot(np.linspace(0, 12, 5), (70,70,70,70,70), dashes=[2, 2], linewidth=3, color='red')
ax.grid()
ax.legend(frameon=False, prop={'size': 14})
plt.savefig('IHE_influence_reinjection_T.pdf')

#%%

fluids = ['R245fa', 'R600', 'R245ca', 'Isopentane']

R245fa = pd.read_csv('IHE_install_R245fa.csv', 
                     header=0, names=['T_1', 'T_2','p_2',
                                      'T_35', 'Q_internal heat exchanger', 'gross power output',
                                      'net power output', 'thermal efficiency', 'net efficiency', 'IHE sizing factor'])
R600 = pd.read_csv('IHE_install_R600.csv', 
                     header=0, names=['T_1', 'T_2','p_2',
                                      'T_35', 'Q_internal heat exchanger', 'gross power output',
                                      'net power output', 'thermal efficiency', 'net efficiency', 'IHE sizing factor'])

fig, ax = plt.subplots(figsize=(8, 6), dpi=100)
ax.plot(R245fa['T_1'], R245fa['net power output'], color='blue', marker="v", label='R245fa')
ax.plot(R600['T_1'], R600['net power output'], color='blue', marker="x", label='R600')
ax.set(xlabel= 'Turbine inlet temperature (°C)', ylabel='Net power output (MW)')
plt.ylim(12.5, 13)
ax2=ax.twinx()
ax2.plot(R245fa['T_1'], R245fa['Q_internal heat exchanger'], color='black', marker="*", label='R245fa')
ax2.plot(R600['T_1'], R600['Q_internal heat exchanger'], color='black', marker="o", label='R600')
ax2.set(ylabel='Required IHE heat exchange rate (MW)')
plt.ylim(5.5, 8)
plt.xlim(125, 131)
ax.yaxis.label.set_color('blue')
ax.yaxis.label.set_size(18)
ax.xaxis.label.set_size(18)
ax.tick_params(axis="x", labelsize=15)
ax.tick_params(axis="y", labelsize=15)
ax.tick_params(axis='y', colors='blue')
ax2.yaxis.label.set_size(17)
ax2.tick_params(axis="y", labelsize=15)
ax.grid()
ax.legend(loc="lower left", frameon=False, prop={'size': 15})
ax2.legend(loc="best", frameon=False, prop={'size': 15})
plt.savefig('IHE_install_R245fa_R600.pdf')


















