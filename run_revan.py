import os
import shutil
import argparse
import yaml

# Initialize parser & read arguments from command line
def parse_args():

	parser = argparse.ArgumentParser()
	parser.add_argument("-y", "--yaml", help = "Path to input .yaml file", default='examples/example_revan.yaml')
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
	mass_model_path = os.path.abspath(inputs['mass_model_path']) + '/'
	config_file_path = os.path.abspath(inputs['config_file_path']) + '/'
	if not os.path.isdir(output_path):
		os.makedirs(output_path)
	if not os.path.isdir(output_path + 'output_files/'):
		os.makedirs(output_path + 'output_files/')

	return input_path, output_path, mass_model_path, config_file_path

# Run revan on all files in input path
def run_revan(input_path, output_path, mass_model_path, config_file_path):

	cwd = os.getcwd()

	for file in os.listdir(input_path):
		if os.path.isfile(input_path + file) and file.split('.')[-2] == 'sim' and not os.path.exists(output_path + file.split('.')[0] + '.inc1.id1.tra.gz'):
			print('Reconstructing ' + file)
			os.chdir(input_path)
			os.system('revan -f ' + file + ' -g ' + mass_model_path + ' -c ' + config_file_path + ' -n -a > ' + output_path + 'output_files/output_' + file.split('.')[0] + '.txt')
			os.chdir(cwd)
			os.system('mv ' + input_path + file.split('.')[0] + '.inc1.id1.tra.gz ' + output_path + file.split('.')[0] + '.inc1.id1.tra.gz')

	if os.path.exists(cwd + '/absorptions'):
		shutil.rmtree(cwd + '/absorptions')


def main():
	input_file = parse_args().yaml
	inputs = read_yaml(input_file)
	input_path, output_path, mass_model_path, config_file_path = define_paths(inputs)
	run_revan(input_path, output_path, mass_model_path, config_file_path)

if __name__ == "__main__":
    main()