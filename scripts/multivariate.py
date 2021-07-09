from geothermal_orc import multivariate_optimization
import json
import sys

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/multivariate.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result, opt_results = multivariate_optimization(**input_data)
