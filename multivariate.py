import matplotlib.pyplot as plt
from geothermal_orc_design import multivariate_optimization
from collections import OrderedDict
from znes_plotting import plot, shared

import itertools

import json
import sys


cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/test_multi.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = multivariate_optimization(**input_data)

for fluid, data in result.items():
    combinations = list(itertools.combinations(input_data['variables'].keys(), 2))

    for combination in combinations:

        ax = plot.scatter(
            data=data[data['valid']], x=combination[0], y=combination[1],
            xlabel=input_data['variables'][combination[0]]['label'],
            ylabel=input_data['variables'][combination[1]]['label'],
            colormap=plt.cm.get_cmap('RdYlBu'), c=input_data['objective'])

    # Create multipage PDF from plots
    shared.create_multipage_pdf(fluid + '_multivariate_scatterplots.pdf')
