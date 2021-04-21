import matplotlib.pyplot as plt
from geothermal_orc_design import single_optimization
from collections import OrderedDict

import itertools

import json
import sys


cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/test.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_optimization(**input_data)

print(result)
