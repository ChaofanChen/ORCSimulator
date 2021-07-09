from geothermal_orc import single_optimization
import json
import sys

cur_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

with open(cur_dir + '/Single_optimization.json', 'r') as f:
    input_data = json.load(f)
    f.close()

result = single_optimization(**input_data)
