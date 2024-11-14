import sys

sys.path.append('/path/to/grb-simulations/') # or save scripts directly in grb_simulator directory

from grb_simulator import parse_args, read_yaml, define_paths, load_megalib, trigger_algorithm_inputs

# Define input yaml file
input_file = parse_args([['-y', '--yaml', 'Path to input .yaml file']]).yaml

# Create trigger algorithm input files
trigger_algorithm = trigger_algorithm_inputs(input_file)
trigger_algorithm.create_event_files()
