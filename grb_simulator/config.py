import os
import yaml
import argparse
import h5py
import csv
import numpy as np
import matplotlib.pyplot as plt

def parse_args(arg_list):
	"""
	Read arguments from command line.

	Parameters
	----------
	arg_list : list of str
		2D list with each row corresponding to an argument and three columns (-flag, --name, help description) 

	Returns
	----------
	args : argparse.Namespace
		Arguments from command line
	"""

	parser = argparse.ArgumentParser()
	for i in range(len(arg_list)):
		parser.add_argument(arg_list[i][0], arg_list[i][1], help=arg_list[i][2])
	args = parser.parse_args()

	return args

def read_yaml(file):
	"""
	Read yaml file.

	Parameters
	----------
	file : str
		Name of yaml file 

	Returns
	----------
	inputs : dict
		Contents of yaml file
	"""

	with open(file, "r") as myfile:
		inputs = yaml.safe_load(myfile)

	return inputs

def fill_dict(keys, data):
	"""
	Fill dictionary with contents.

	Parameters
	----------
	keys : list of str
		Dictionary keys 
	data : np.array
		Data to add to dictionary

	Returns
	----------
	this_dict : dict
		Output dictionary
	"""

	this_dict = {}
	for i in range(len(keys)):
		if isinstance(data[i], float):
			this_dict[keys[i]] = float(data[i])
		else:
			try:
				this_dict[keys[i]] = int(data[i])
			except:
				this_dict[keys[i]] = str(data[i])

	return this_dict
	
def write_yaml(file, yaml_dict):
	"""
	Write yaml file.

	Parameters
	----------
	file : str
		Name of yaml file 
	yaml_dict : dict
		Contents of yaml file 
	"""

	with open(file, 'w') as f:
		yaml.dump(yaml_dict, f, sort_keys=False)

def define_paths(paths, create_dir):
	"""
	Define paths.

	Parameters
	----------
	paths : list of str
		List of paths
	create_dir : list of bool
		Whether to create directory for each path

	Returns
	----------
	output_paths : list of str
		List of paths with proper formatting
	"""

	output_paths = []

	if len(paths) == len(create_dir):
		for i in range(len(paths)):
			this_path = paths[i]
			if create_dir[i]:
				if not os.path.isdir(this_path):
					os.makedirs(this_path)
			if os.path.isfile(this_path):
				this_path = os.path.abspath(paths[i])
			elif os.path.isdir(this_path):
				this_path = os.path.abspath(paths[i]) + '/'
			else:
				raise RuntimeError(this_path + ' does not exist')
			output_paths.append(this_path)

	return output_paths

def read_hdf5(file):
	"""
	Read trigger algorithm hdf5 file.

	Parameters
	----------
	file : str
		Name of file 
	"""

	with h5py.File(file, 'r') as hf:
		time_array = hf['trigger_data'][0]
		energy_array = hf['trigger_data'][1]

	return time_array, energy_array

def read_event_list(event_list):
	"""
	Read event list file and save as a dictionary.

	Parameters
	----------
	event_list : str
		Event list data

	Returns
	----------
	data : dict
		Dictionary containing event list
	"""

	with open(event_list, 'r', newline='', encoding='utf-8') as csvfile:
		reader = csv.reader(csvfile)
		header = next(reader)  # Read the header row
		data = {}
		for col_name in header:
			data[col_name] = []

		for row in reader:
			for i, value in enumerate(row):
				data[header[i]].append(value)

	return data

def plot_triggers(path, plot_path, bin_size=0.05, time_range=(1835487300., 1835573700.)):
	"""
	Plot lightcurve for each event.

	Parameters
	----------
	path : str
		Path to trigger algorithm input directory
	plot_path : str
		Path to directory to save plots
	bin_size : float, optional
		Time bin size in s
	time_range : tuple of floats, optional
		Start and end times  in s
	"""

	nbins = int((time_range[1] - time_range[0]) / bin_size)
	lightcurves = {}

	colors = ['red', 'green', 'blue', 'orange', 'purple', 'pink']

	for file in os.listdir(path):
		print(file)
		if file.endswith('hdf5'):
			time_array, energy_array = read_hdf5(path + file)
			mask = (energy_array >= 80.) & (energy_array <= 2000.)
			time_array_filtered = time_array[mask]
			counts, bins = np.histogram(time_array_filtered, bins=nbins, range=time_range)
			lightcurves[file[:5]] = counts
		elif file.endswith('txt'):
			with open(path + file, 'r', newline='', encoding='utf-8') as csvfile:
				reader = csv.reader(csvfile, delimiter='\t')
				header = next(reader)
				data = {}
				for col_name in header:
					data[col_name] = []
				for row in reader:
					for i, value in enumerate(row):
						data[header[i]].append(value)

	for i in range(len(data['Event name'])):
		print(data['Event name'][i])
		start_time = float(data['Start time (s)'][i])
		end_time = start_time + float(data['Duration (s)'][i])
		plot_range = [start_time - 30, start_time + 40]
		start_index = np.argmin(np.abs(bins - plot_range[0]))
		end_index = np.argmin(np.abs(bins - plot_range[1]))
		these_bins = bins[start_index:end_index]
		for j, key in enumerate(lightcurves.keys()):
			these_counts = lightcurves[key][start_index:end_index - 1]
			plt.stairs(these_counts, these_bins, color=colors[j], alpha=0.4)
		plt.title(data['Event name'][i] + ' (T90 = ' + data['GBM T90 (s)'][i] + ' s)')
		plt.axvspan(start_time, end_time, facecolor='gray', alpha=0.5)
		plt.savefig(plot_path + data['Event name'][i] + '.png', dpi=350)
		plt.clf()

class suppress_output(object):
    
  def __init__(self):
    """
    Suppress command line output.
    """

    self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)] # Open a pair of null files
    self.save_fds = [os.dup(1), os.dup(2)] # Save actual stdout (1) and stderr (2) file descriptors

  def __enter__(self):
    """
    Assign null pointers to stdout and stderr.
    """

    os.dup2(self.null_fds[0],1)
    os.dup2(self.null_fds[1],2)

  def __exit__(self, *_):
    """
    Re-assign real stdout & stderr back to (1) & (2).
    """

    os.dup2(self.save_fds[0],1)
    os.dup2(self.save_fds[1],2)

    for fd in self.null_fds + self.save_fds: # Close all file descriptors
      os.close(fd)
