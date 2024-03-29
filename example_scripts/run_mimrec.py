import sys

sys.path.append('/path/to/grb-simulations/') # or save scripts directly in grb_simulator directory

from grb_simulator import parse_args, read_yaml, run_megalib

# Read input yaml file
input_file = parse_args([['-y', '--yaml', 'Path to input .yaml file']]).yaml
inputs = read_yaml(input_file)

# Run mimrec on all files in input path
megalib = run_megalib(inputs)
megalib.run_mimrec()