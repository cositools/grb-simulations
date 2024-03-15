import sys

sys.path.append('/path/to/grb-simulations/') # or save scripts directly in grb_simulator directory

from grb_simulator import parse_args, read_yaml, define_paths, source_files

# Read input yaml file
input_file = parse_args([['-y', '--yaml', 'Path to input .yaml file']]).yaml
inputs = read_yaml(input_file)

# Create MEGAlib source files for all events in input directory
create_files = source_files(inputs)
create_files.write_readme()
create_files.make_source_files()