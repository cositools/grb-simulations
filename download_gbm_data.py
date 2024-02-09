from gbm.finder import BurstCatalog
from gbm.finder import TriggerFtp
import os
import yaml
import argparse
import numpy as np

# Initialize parser & read arguments from command line
def parse_args():

	parser = argparse.ArgumentParser()
	parser.add_argument("-y", "--yaml", help = "Path to input .yaml file", default='examples/example_download.yaml')
	args = parser.parse_args()

	return args

# Read yaml file
def read_yaml(file):

	with open(file, "r") as myfile:
		inputs = yaml.safe_load(myfile)

	return inputs

# Define path to output files
def define_paths(inputs):

	output_path = os.path.abspath(inputs['output_path']) + '/'
	if not os.path.isdir(output_path):
		os.makedirs(output_path)

	return output_path

# Define parameters used to filter bursts
def define_parameters(inputs, burst_catalog):

	filter_list = []
	try:
		for item in inputs['filters']:
			filter_list.append((item, inputs['filters'][item][0], inputs['filters'][item][1]))
		sliced_burst_catalog = burst_catalog.slices(filter_list)
	except:
		print('No filters defined. Downloading all GBM trggers')
		sliced_burst_catalog = burst_catalog


	return sliced_burst_catalog

# Write yaml file with data for each burst
def write_yaml(inputs, data, output_path):

	mydict = {}
	for i in range(len(inputs['download'])):
		if type(data[i]) == np.float64:
			mydict[inputs['download'][i]] = float(data[i])
		else:
			try:
				mydict[inputs['download'][i]] = int(data[i])
			except:
				mydict[inputs['download'][i]] = str(data[i])
	with open(output_path + data[inputs['download'].index('trigger_name')] + '/' + data[inputs['download'].index('trigger_name')] + '.yaml', 'w') as f:
		yaml.dump(mydict, f)

# Download tte data
def download_bursts(inputs, burst_catalog, output_path):

	sliced_burst_catalog = define_parameters(inputs, burst_catalog)
	burst_list = sliced_burst_catalog.get_table(columns=inputs['download'])
	trig_finder = TriggerFtp(burst_list[len(burst_list)-1][inputs['download'].index('trigger_name')].replace('bn', ''))

	for i in range(len(burst_list)):
		name = burst_list[i][inputs['download'].index('trigger_name')]
		print('Downloading data for ' + name + ' (' + str(i+1) + '/' + str(len(burst_list)) + ')')
		if not os.path.isdir(output_path + name):
			os.mkdir(output_path + name)
		write_yaml(inputs, burst_list[i], output_path)
		detectors = burst_list[i][inputs['download'].index('bcat_detector_mask')]
		trig_finder.set_trigger(name.replace('bn', ''))
		
		detector_list = []
		for j in range(len(detectors)):
			if detectors[j] == '1':
				if j <= 9:
					det = 'n' + str(j)
				elif j == 10:
					det = 'na'
				elif j == 11:
					det = 'nb'
				detector_list.append(det)
		bgo1 = 0
		bgo2 = 0
		for item in detectors[0:6:1]:
			bgo1 +=float(item)
		for item in detectors[6:13:1]:
			bgo2 +=float(item)
		if bgo1 > bgo2:
			detector_list.append('b0')
		elif bgo1 < bgo2:
			detector_list.append('b1')
		else:
			detector_list.append('b0')
			detector_list.append('b1')
		trig_finder.get_tte(output_path + name, dets=detector_list)


def main():
	input_file = parse_args().yaml
	inputs = read_yaml(input_file)
	output_path = define_paths(inputs)
	burst_catalog = BurstCatalog()
	sliced_burst_catalog = download_bursts(inputs, burst_catalog, output_path)

if __name__ == "__main__":
    main()
