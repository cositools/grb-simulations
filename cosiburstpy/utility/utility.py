import os
import yaml
import argparse
import h5py
import csv
from datetime import date
import logging
from cosiburstpy._version import version

logger = logging.getLogger(__name__)

def parse_args(arg_list):
	'''
	Read arguments from command line.

	Parameters
	----------
	arg_list : list of str
		2D list with each row corresponding to an argument and three columns (-flag, --name, help description) 

	Returns
	-------
	args : argparse.Namespace
		Arguments from command line
	'''

	parser = argparse.ArgumentParser()

	for i in range(len(arg_list)):
		parser.add_argument(arg_list[i][0], arg_list[i][1], help=arg_list[i][2])

	args = parser.parse_args()

	return args

def read_yaml(file):
	'''
	Read .yaml file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to .yaml file 

	Returns
	-------
	data : dict
		Contents of .yaml file
	'''

	logger.info(f'Reading file: {file}')

	with open(file, 'r') as f:
		data = yaml.safe_load(f)

	return data

def write_yaml(file, data):
	'''
	Write .yaml file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to .yaml file 
	data : dict
		Contents of .yaml file
	'''

	logger.info(f'Writing file: {file}')

	with open(file, 'w') as f:
		yaml.dump(data, f, sort_keys=False)

def read_hdf5(file):
	'''
	Read .hdf5 file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to .hdf5 file

	Returns
	-------
	data : dict
		Contents of .hdf5 file
	file_attributes : dict
		Attributes of .hdf5 file
	dataset_attributes : dict
		Attributes of datasets in .hdf5 file
	'''

	logger.info(f'Reading file: {file}')

	data = {}
	file_attributes = {}
	dataset_attributes = {}

	with h5py.File(file, 'r') as f:

		for key in f.attrs:
			file_attributes[key] = f.attrs[key]

		for key in f:

			dataset = []

			for item in f[key]:
				dataset.append(item)

			data[key] = dataset

			for att_key in f[key].attrs:
				dataset_attributes[key][att_key] = f[key].attrs[att_key]

	if file_attributes == {}:
		file_attributes = None

	if dataset_attributes == {}:
		dataset_attributes = None

	return data, file_attributes, dataset_attributes

def write_hdf5(file, data, file_attributes=None, dataset_attributes=None):
	'''
	Write .hdf5 file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to .hdf5 file
	data : dict
		Contents of .hdf5 file
	file_attributes : dict
		Attributes of .hdf5 file
	dataset_attributes : dict
		Attributes of datasets in .hdf5 file
	'''

	logger.info(f'Writing file: {file}')

	with h5py.File(file, 'w') as f:

		if file_attributes:

			for key in file_attributes:
				f.attrs[key] = file_attributes[key]

		for key in data:

			dataset = f.create_dataset(key, data=data[key], compression='gzip')

			if dataset_attributes and key in dataset_attributes:
				for att_key in dataset_attributes[key]:
					dataset.attrs[key][att_key] = dataset_attributes[key][att_key]

def write_readme(file, inputs_path=None, input_parameters=None):
	'''
	Write README file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to README file
	inputs_path : str
		Path to directory containing script and configuration file run to produce output
	input_parameters : dict
		Input parameters
	'''

	logger.info(f'Writing file: {file}')

	with open(file, 'w') as f:

		f.write(f'This output was generated on {date.today()} using cosipyburst version {version}')

		if inputs_path:
			f.write(f' by running the files located in {inputs_path}.  \n')

		else:
			f.write(f'.  ')

		if input_parameters:

			f.write(f'\n\nParameters used:  \n')

			for key in input_parameters:
				f.write(f'{key}: {input_parameters[key]}  \n')

def read_csv(file, delimiter=None, ignore_spaces=True, encoding='utf-8'):
	'''
	Read .csv file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to .csv file 
	delimiter : str, optional
		Delimiter for .csv file
	skipinitialspace ; bool, optional
		Whether to ignore spaces after delimiter
	encoding : str, optional
		Encoding of .csv file

	Returns
	-------
	data : dict
		Contents of .csv file
	'''

	logger.info(f'Reading file: {file}')

	data = {}

	with open(file, 'r', newline='', encoding=encoding) as f:

		if delimiter:
			reader = csv.reader(f, delimiter=delimiter, skipinitialspace=ignore_spaces)

		else:
			reader = csv.reader(f, skipinitialspace=ignore_spaces)

		header = next(reader)

		for column in header:
			data[column] = []

		for row in reader:
			for i, value in enumerate(row):
				data[header[i]].append(value)

	return data

def write_csv(file, data, delimiter='\t'):
	'''
	Write .csv file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to .csv file
	data : dict
		Contents of .csv file
	delimiter : str, optional
		Delimiter for .csv file
	'''

	logger.info(f'Writing file: {file}')

	with open(file, 'w', newline='') as f:

		w = csv.writer(f, delimiter=delimiter)
		w.writerow(data.keys())
		w.writerows(zip(*data.values()))

def make_dict(keys, contents):
	'''
	Fill dictionary with contents.

	Parameters
	----------
	keys : list of str
		Dictionary keys 
	contents : list
		Contents to add to dictionary

	Returns
	-------
	data : dict
		Data dictionary
	'''

	data = {}

	for i, key in enumerate(keys):

		try:

			value = float(contents[i])

			if int(value) == value:
				data[key] = int(value)

			else:
				data[key] = value

		except:

			data[key] = str(contents[i])

	return data

def search_and_replace(input_file, output_file, search, replace):
	'''
	Search and replace within text file.

	Parameters
	----------
	input_file : pathlib.PosixPath
		Path to original file
	output_file : pathlib.PosixPath
		Path to updated file
	search : str
		Text to be replaced
	replace : str
		Text to replace
	'''

	with open(input_file, 'r') as f:
		input_data = f.read()
	
	output_data = input_data.replace(search, replace)

	with open(output_file, 'w') as f:
		f.write(output_data)

class SuppressOutput(object):
	
	def __init__(self):
		'''
		Suppress command line output.
		'''

		self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)] # Open a pair of null files
		self.save_fds = [os.dup(1), os.dup(2)] # Save actual stdout (1) and stderr (2) file descriptors

	def __enter__(self):
		'''
		Assign null pointers to stdout and stderr.
		'''

		os.dup2(self.null_fds[0],1)
		os.dup2(self.null_fds[1],2)

	def __exit__(self, *_):
		'''
		Re-assign real stdout & stderr back to (1) & (2).
		'''

		os.dup2(self.save_fds[0],1)
		os.dup2(self.save_fds[1],2)

		# Close all file descriptors
		for fd in self.null_fds + self.save_fds:
			os.close(fd)
