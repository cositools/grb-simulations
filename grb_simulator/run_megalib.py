import os
import shutil
import glob
import fnmatch
from .config import read_yaml, define_paths

class run_megalib():

	def __init__(self, input_file={}):
		"""
		Run cosima, revan, and/or mimrec.

		Parameters
		----------
		input_file : str
			Path to input .yaml file
		"""

		inputs = read_yaml(input_file)


		if 'source_files' in inputs['paths'] and 'sim_files' in inputs['paths']:
			[self.input_cosima, self.output_cosima] = define_paths([inputs['paths']['source_files'], inputs['paths']['sim_files']], [False, True])

			if 'parallel' in inputs['cosima'] and inputs['cosima']['parallel']:
				self.cosima_parallel = True
				if 'instances' in inputs['cosima']:
					self.cosima_instances = int(inputs['cosima']['instances'])
			else:
				self.cosima_parallel = False

			if 'zip' in inputs['cosima'] and not inputs['cosima']['zip']:
				self.cosima_zip = False
			else:
				self.cosima_zip = True

		if 'sim_files' in inputs['paths'] and 'tra_files' in inputs['paths'] and 'mass_model' in inputs['paths'] and 'config' in inputs['revan']:
			[self.input_revan, 
			 self.output_revan, 
			 mass_model_path,
			 self.config_revan] = define_paths([inputs['paths']['sim_files'], inputs['paths']['tra_files'], 
			 									inputs['paths']['mass_model'], inputs['revan']['config']], 
			 								   [False, True, False, False])

		if 'tra_files' in inputs['paths'] and 'extracted_tra_files' in inputs['paths'] and 'mass_model' in inputs['paths'] and 'config' in inputs['mimrec']:
			[self.input_mimrec, 
			 self.output_mimrec, 
			 mass_model_path,
			 self.config_mimrec] = define_paths([inputs['paths']['tra_files'], inputs['paths']['extracted_tra_files'], 
			 									 inputs['paths']['mass_model'], inputs['mimrec']['config']], 
			 								    [False, True, False, False])

	def cosima(self, file, output_path, zipped):
		"""
		Run cosima on a source file.

		Parameters
		----------
		file : str
			Source file name 
		"""

		if not os.path.isdir(output_path + 'output/'):
			os.makedirs(output_path + 'output/')

		if zipped:
			os.system('cosima ' + file + ' > ' + output_path + 'output/output_' + file.split('.')[0] + '.txt')
		else:
			os.system('cosima -u ' + file + ' > ' + output_path + 'output/output_' + file.split('.')[0] + '.txt')

	def mcosima(self, file, output_path, zipped, instances=None):
		"""
		Run mcosima on a source file.

		Parameters
		----------
		file : str
			Source file name 
		"""

		if not os.path.isdir(output_path + 'output/'):
			os.makedirs(output_path + 'output/')

		if instances == None:
			if zipped:
				os.system('mcosima -z ' + file + ' > ' + output_path + 'output/output_' + file.split('.')[0] + '.txt')
			else:
				os.system('mcosima ' + file + ' > ' + output_path + 'output/output_' + file.split('.')[0] + '.txt')
		else:
			if zipped:
				os.system('mcosima -t ' + instances + ' -z ' + file + ' > ' + output_path + '/output_' + file.split('.')[0] + '.txt')
			else:
				os.system('mcosima -t ' + instances + file + ' > ' + output_path + '/output_' + file.split('.')[0] + '.txt')

	def revan(self, file, output_path, mass_model, config_file):
		"""
		Run revan on a sim file.

		Parameters
		----------
		file : str
			Sim or sim.gz file name 
		"""

		if not os.path.isdir(output_path + 'output/'):
			os.makedirs(output_path + 'output/')
		os.system('revan -f ' + file + ' -g ' + mass_model + ' -c ' + config_file + ' -n -a > ' + output_path + 'output/output_' + file.split('.')[0] + '.txt')

	def mimrec(self, file, output_path, mass_model, config_file):
		"""
		Run mimrec on a tra file.

		Parameters
		----------
		file : str
			Tra or tra.gz file name 
		"""

		if not os.path.isdir(output_path + 'output/'):
			os.makedirs(output_path + 'output/')
		os.system('mimrec -f ' + file + ' -g ' + mass_model + ' -c ' + config_file + ' -n -x > ' + output_path + 'output/output_' + file.split('.')[0] + '.txt')

	def run_on_all(self, run_func, input_path, output_path, input_ext, output_ext, print_text):
		"""
		Run cosima, mcosima, revan, or mimrec on all files in a directory.

		Parameters
		----------
		run_func : function
			Function to run (run_cosima, run_mcosima, run_revan, or run_mimrec)
		input_path : str
			Path to directory housing input files
		output_path : str
			Path to directory housing output files
		input_ext : list of str
			File extensions of input files
		output_ext : str
			File extension of output files
		print_text : str
			Text to print when running
		"""

		file_list = glob.glob(output_path + '*.' + output_ext)
		for i in range(len(file_list)):
			file_list[i] = file_list[i].split('/')[-1]

		cwd = os.getcwd()
		for file in os.listdir(input_path):
			for ext in input_ext:
				if os.path.isfile(input_path + file) and file[-len(ext):] == ext:
					if not (len(fnmatch.filter(file_list, file.split('.')[0] + '*' + output_ext)) > 0):
						print(print_text + ' ' + file)
						os.chdir(input_path)
						run_func(file)
						os.chdir(cwd)
						os.system('mv ' + input_path + '*.' + output_ext + ' ' + output_path)
					else:
						print("Output file for " + file + " already exists.")

		if os.path.exists(cwd + '/absorptions'):
			shutil.rmtree(cwd + '/absorptions')

	def run_cosima(self, input_path=None, output_path=None, zipped=None, instances=None):
		"""
		Run cosima or mcosima on all source files in a directory.

		Parameters
		----------
		input_path : str, optional
			Path to directory housing source files, if different from path defined when creating run_megalib object
		output_path : str, optional
			Path to directory to store output sim or sim.gz files, if different from path defined when creating run_megalib object
		zipped : bool, optional
			Whether to zip output .sim files
		instances : int, optional
			Number of instances to run for mcosima, default is number of CPU cores
		"""

		if input_path == None:
			if hasattr(self, 'input_cosima'):
				input_path = self.input_cosima
			else:
				raise RuntimeError(".source file directory must be provided in input .yaml file or as an argument when running run_cosima().")
		if output_path == None:
			if hasattr(self, 'output_cosima'):
				output_path = self.output_cosima
			else:
				raise RuntimeError(".sim file directory must be provided in input .yaml file or as an argument when running run_cosima().")
		if zipped == None:
			zipped = self.cosima_zip

		if not self.cosima_parallel:
			cosima_func = lambda file: self.cosima(file, output_path, self.cosima_zip)
		else:
			if instances == None and hasattr(self, 'cosima_instances'):
				instances = self.cosima_instances
			cosima_func = lambda file: self.mcosima(file, output_path, self.cosima_zip, instances=instances)
			
		if zipped:
			self.run_on_all(mcosima_func, input_path, output_path, ['source'], 'sim.gz', 'Simulating')
		else:
			self.run_on_all(mcosima_func, input_path, output_path, ['source'], 'sim', 'Simulating')

		with open(output_path + 'README.md', 'w') as f:
			f.write('# .sim Files\n\n')
			f.write('This directory contains .sim files created using cosima from the .source files in `' + input_path + '`.')
			if self.cosima_parallel:
				if instances == None:
					f.write('This directory contains .sim files created using mcosima from the .source files in `' + input_path + '`.')
				else:
					f.write('This directory contains .sim files created using mcosima with ' + str(instances) + ' instances from the .source files in `' + input_path + '`.')

	def run_revan(self, input_path=None, output_path=None, config_file=None):
		"""
		Run revan on all sim and sim.gz files in a directory.

		Parameters
		----------
		input_path : str, optional
			Path to directory housing sim and/or sim.gz files, if different from path defined when creating run_megalib object
		output_path : str, optional
			Path to directory to store output tra or tra.gz files, if different from path defined when creating run_megalib object
		config_file : str, optional
			Path to revan configuration file, if different from path defined when creating run_megalib object
		zipped : bool, optional
			Whether to zip output tra files
		"""

		if input_path == None:
			if hasattr(self, 'input_revan'):
				input_path = self.input_revan
			else:
				raise RuntimeError(".sim file directory must be provided in input .yaml file or as an argument when running run_revan().")
		if output_path == None:
			if hasattr(self, 'output_revan'):
				output_path = self.output_revan
			else:
				raise RuntimeError(".tra file directory must be provided in input .yaml file or as an argument when running run_revan().")
		if config_file == None:
			if hasattr(self, 'config_revan'):
				config_file = self.config_revan
			else:
				raise RuntimeError("Configuration file path must be provided in input .yaml file or as an argument when running run_revan().")
		if not hasattr(self, 'mass_model'):
			raise RuntimeError("Mass model path must be provided in input .yaml file or as an argument when running run_revan().")

		revan_func = lambda file: self.revan(file, output_path, self.mass_model, config_file)

		self.run_on_all(revan_func, input_path, output_path, ['sim','sim.gz'], 'tra.gz', 'Reconstructing')

		with open(output_path + 'README.md', 'w') as f:
			f.write('# .tra Files\n\n')
			f.write('This directory contains .tra files created using revan from the .sim files in `' + input_path + '` using `' + config_file + '` as the configuration file.')

	def run_mimrec(self, input_path=None, output_path=None, config_file=None):
		"""
		Run mimrec on all tra and tra.gz files in a directory.

		Parameters
		----------
		input_path : str, optional
			Path to directory housing tra and/or tra.gz files, if different from path defined when creating run_megalib object
		output_path : str, optional
			Path to directory to store output extracted tra files, if different from path defined when creating run_megalib object
		config_file : str, optional
			Path to mimrec configuration file, if different from path defined when creating run_megalib object
		"""

		if input_path == None:
			if hasattr(self, 'input_mimrec'):
				input_path = self.input_mimrec
			else:
				raise RuntimeError(".tra file directory must be provided in input .yaml file or as an argument when running run_mimrec().")
		if output_path == None:
			if hasattr(self, 'output_mimrec'):
				output_path = self.output_mimrec
			else:
				raise RuntimeError("Extracted .tra file directory must be provided in input .yaml file or as an argument when running run_mimrec().")
		if config_file == None:
			if hasattr(self, 'config_mimrec'):
				config_file = self.config_mimrec
			else:
				raise RuntimeError("Configuration file path must be provided in input .yaml file or as an argument when running run_mimrec().")
		if not hasattr(self, 'mass_model'):
			raise RuntimeError("Mass model path must be provided in input .yaml file or as an argument when running run_mimrec().")

		mimrec_func = lambda file: self.mimrec(file, output_path, self.mass_model, config_file)

		self.run_on_all(mimrec_func, input_path, output_path, ['tra','tra.gz'], 'extracted.tra', 'Extracting events for')

		with open(output_path + 'README.md', 'w') as f:
			f.write('# Extracted .tra Files\n\n')
			f.write('This directory contains extracted .tra files created using mimrec from the .tra files in `' + input_path + '` using `' + config_file + '` as the configuration file.')