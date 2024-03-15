import sys

sys.path.append('/path/to/grb-simulations/') # or save scripts directly in grb_simulator directory

from grb_simulator import parse_args, read_yaml, download_gbm_data

# Read input yaml file
input_file = parse_args([['-y', '--yaml', 'Path to input .yaml file']]).yaml
inputs = read_yaml(input_file)

# Download data from GBM
gbm_download = download_gbm_data(inputs)
gbm_download.download_events()