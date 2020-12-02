#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 15:56:49 2020

@author: chencha
"""

import matplotlib.pyplot as plt
import geothermal_orc_design_fwi_chaofan 

# -----------parametric optimization for every working fluid-------------------------------------------
import pygmo as pg

class optimization_problem():

    def fitness(self, x):
        f1 = 1 / self.model.calculate_efficiency_opt_with_ihe(x, working_fluid='R600')
        ci1 = -x[0] + x[1]
        print(x)
        return [f1, ci1]

    def get_nobj(self):
        """Return number of objectives."""
        return 1

    # equality constraints
    def get_nec(self):
        return 0

    # inequality constraints
    def get_nic(self):
        return 1

    def get_bounds(self):
        """Return bounds of decision variables."""
        return ([12, -10], [16, -7])

optimize = optimization_problem()
optimize.model = geothermal_orc_design_fwi_chaofan.PowerPlant(working_fluid='R600')
prob = pg.problem(optimize)
num_gen = 15

pop = pg.population(prob, size=10)
algo = pg.algorithm(pg.ihs(gen=num_gen))

result = {'champion': [], 'power_output': [], 'generation': [],
          'Turbine inlet pressure': [], 'Q_ihe': []}

for gen in range(num_gen):
    result["generation"].append(gen)
    result["champion"].append(1 / pop.champion_f[0])

    decision_var = pop.get_x()
    for Td_bp in decision_var:
        result['Turbine inlet pressure'].append(Td_bp[0])
        result['Q_ihe'].append(Td_bp[1])

    fitness = pop.get_f()
    for power_output in fitness:
        result['power_output'].append(1/power_output[0])

    print()
    print('Evolution: {}'.format(gen))
    print('Power Output: {} MW'.format(round(1 / pop.champion_f[0], 4)))
    pop = algo.evolve(pop)

print()
print('Power Output: {} MW'.format(round(1 / pop.champion_f[0], 4)))
print('Turbine inlet pressure: {} bar'.format(round(pop.champion_x[0], 4)))
print('Q_ihe: {} W'.format(round(pop.champion_x[1], 4)))

# scatter plot
cm = plt.cm.get_cmap('RdYlBu')
sc = plt.scatter(result['Q_ihe'], result['Turbine inlet pressure'], linewidth=0.25,
                 c=result['power_output'], cmap=cm, alpha=0.5, edgecolors='black')
plt.scatter(pop.champion_x[1], pop.champion_x[0], marker='x', linewidth=1,
            c='red')
plt.annotate('Optimum', xy=(pop.champion_x[1], pop.champion_x[0]),
             xytext=(pop.champion_x[1]+0.0002, pop.champion_x[0]+0.0002),
             arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',
                             color='red')
             )
plt.ylabel('Turbine inlet pressure [bar]')
plt.xlabel('Q_ihe [MW]')
plt.colorbar(sc, label='Power Output [MW]')
plt.savefig("P_tur_Q_ihe_optimization_power_output.png")
plt.show()