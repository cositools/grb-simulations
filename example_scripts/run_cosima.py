import sys

sys.path.append('/path/to/grb-simulations/') # or save scripts directly in grb_simulator directory

from grb_simulator import parse_args, read_yaml, run_megalib

# Define input yaml file
input_file = parse_args([['-y', '--yaml', 'Path to input .yaml file']]).yaml
inputs = read_yaml(input_file)

# Run cosima or mcosima on all files in input path
megalib = run_megalib(input_file)
if 'parallel' in inputs.keys() and inputs['parallel']:
	if 'instances' in inputs.keys():
		megalib.run_mcosima(instances=str(inputs['instances']))
	else:
		megalib.run_mcosima()
else:	
	megalib.run_cosima()