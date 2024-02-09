import os
import shutil
import argparse
import yaml

# Initialize parser & read arguments from command line
def parse_args():

	parser = argparse.ArgumentParser()
	parser.add_argument("-y", "--yaml", help = "Path to input .yaml file", default='examples/example_cosima.yaml')
	args = parser.parse_args()

	return args

# Read yaml file
def read_yaml(file):

	with open(file, "r") as myfile:
		inputs = yaml.safe_load(myfile)

	return inputs

# Define paths to input files and output source files
def define_paths(inputs):

	input_path = os.path.abspath(inputs['input_path']) + '/'
	output_path = os.path.abspath(inputs['output_path']) + '/'
	if not os.path.isdir(output_path):
		os.makedirs(output_path)
	if not os.path.isdir(output_path + 'output_files/'):
		os.makedirs(output_path + 'output_files/')

	return input_path, output_path

# Run cosima on all files in input path
def run_cosima(input_path, output_path):

	cwd = os.getcwd()

	for file in os.listdir(input_path):
		if file.split('.')[-1] == 'source' and not (os.path.exists(output_path + file.split('.')[0] + '.inc1.id1.sim.gz')):
			print('Simulating ' + file)
			os.chdir(input_path)
			os.system('cosima -z ' + file + ' > ' + output_path + 'output_files/output_' + file.split('.')[0] + '.txt')
			os.chdir(cwd)
			os.system('mv ' + input_path + file.split('.')[0] + '.inc1.id1.sim.gz ' + output_path + file.split('.')[0] + '.inc1.id1.sim.gz')

	if os.path.exists(cwd + '/absorptions'):
		shutil.rmtree(cwd + '/absorptions')


def main():
	input_file = parse_args().yaml
	inputs = read_yaml(input_file)
	input_path, output_path = define_paths(inputs)
	run_cosima(input_path, output_path)

if __name__ == "__main__":
    main()