import sys

sys.path.append('/path/to/grb-simulations/') # or save scripts directly in grb_simulator directory

from grb_simulator import parse_args, read_yaml, gbm_to_megalib_inputs

# Read input yaml file
input_file = parse_args([['-y', '--yaml', 'Path to input .yaml file']]).yaml
inputs = read_yaml(input_file)

# Create MEGAlib input files
events = gbm_to_megalib_inputs(inputs)
events.make_spectra_lightcurves()