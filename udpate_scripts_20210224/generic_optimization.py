#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 13:54:02 2020

@author: chencha
"""
import matplotlib.pyplot as plt
import geothermal_orc_design
#from tespy.tools.helpers import UserDefinedEquation
#from tespy.tools.fluid_properties import h_mix_pT, T_mix_ph
from collections import OrderedDict
#import CoolProp as CP
#import numpy as np
import pandas as pd
# -----------parametric optimization for every working fluid-------------------
import pygmo as pg

from znes_plotting import plot, shared

import itertools


class optimization_problem():

    def fitness(self, x):
        inputs = {'IHE_sizing': 0.999}
        i = 0
        for val in x:
            inputs.update({list(self.params_to_opt.keys())[i]: x[i]})
            i += 1

        self.model.run_simulation(**inputs)
        f1 = self.objective_value()
        return [f1]

    def get_nobj(self):
        """Return number of objectives."""
        return 1

    # equality constraints
    def get_nec(self):
        return 0

    # inequality constraints
    def get_nic(self):
        return 0

    def get_bounds(self):
        """Return bounds of decision variables."""
        return (
            [self.params_to_opt[param]['min'] for param in self.params_to_opt.keys()],
            [self.params_to_opt[param]['max'] for param in self.params_to_opt.keys()])


optimize = optimization_problem()
optimize.model = geothermal_orc_design.PowerPlant(working_fluid='R245CA')
optimize.model.nw.get_comp('internal heat exchanger').set_attr(pr1=.98, pr2=.98)


optimize.params_to_opt = OrderedDict(
    Q_brine_ev={'min': -3e7, 'max': -1, 'unit': 'W', 'label': 'Heat exchange rate of brine evaporator (W)'},
#    IHE_sizing={'min': 0.999, 'max': 0.999, 'unit': '1', 'label': 'Internal heat exchanger size'},
    T_air_hot={'min': 13, 'max': 22, 'unit': 'K', 'label': 'Cooling air temperature difference in condenser (°C)'}
)

optimize.objective_value = optimize.model.get_net_power
objectives_list = ['Net power output (MW)']
params_list = list(optimize.params_to_opt.keys())

result = {'champion': [], 'generation': []}

# run optimization problem
prob = pg.problem(optimize)
num_gen = 30
num_ind = 20
pop = pg.population(prob, size=num_ind)
algo = pg.algorithm(pg.de(gen=1))

individuals = pd.DataFrame(
    columns=params_list + objectives_list,
    index=[(gen, ind) for gen in range(num_gen) for ind in range(num_ind)])

for gen in range(num_gen):
    result["generation"].append(gen)
    result["champion"].append(pop.champion_f[0])

    individual = 0
    for x in pop.get_x():
        for i in range(len(x)):
            individuals.loc[[(gen, individual)], params_list[i]] = x[i]
        individual += 1

    individual = 0
    for objective in pop.get_f():
        for i in range(len(objective)):
            individuals.loc[[(gen, individual)], objectives_list[i]] = -objective[i]/1e6
        individual += 1

    print()
    print('Evolution: {}'.format(gen))
    print('Objective value: {}'.format(round(pop.champion_f[0], 4)))
    pop = algo.evolve(pop)

print()
for i in range(len(objectives_list)):
    print(objectives_list[i] + ': {}'.format(round(pop.champion_f[i], 4)))
for i in range(len(params_list)):
    print(params_list[i] + ': {}'.format(round(pop.champion_x[i], 4)))

combinations = list(itertools.combinations(params_list, 2))

for combination in combinations:

    ax = plot.scatter(
        data=individuals, x=combination[0], y=combination[1],
        xlabel=optimize.params_to_opt[combination[0]]['label'],
        ylabel=optimize.params_to_opt[combination[1]]['label'],
        colormap=plt.cm.get_cmap('RdYlBu_r'), c=objectives_list[0])

# Create multipage PDF from plots
shared.create_multipage_pdf('R245CA_multi_opt.pdf')
