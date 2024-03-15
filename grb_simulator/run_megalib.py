import os
import shutil
import glob
import fnmatch
from .config import define_paths

class run_megalib():

	def __init__(self, inputs={}, input_path=None, output_path=None, mass_model_path=None, config_revan=None, config_mimrec=None):
		"""
		Run cosima, revan, and/or mimrec.

		Parameters
		----------
		inputs : dict, optional
			Contents of input yaml file
		input_path : str, optional
			Path to directory housing MEGAlib input files
		output_path : str, optional
			Path to directory to store MEGAlib output files
		mass_model_path : str, optional
			Path to instrument mass model, required if running revan or mimrec
		config_revan : str, optional
			Path to revan configuration file
		config_mimrec : str, optional
			Path to mimrec configuration file
		"""

		if not input_path == None:
			self.input_cosima = input_path
			self.input_revan = input_path
			self.input_mimrec = input_path

		if not output_path == None:
			self.output_cosima = output_path
			self.output_revan = output_path
			self.output_mimrec = output_path

		if not mass_model_path == None:
			self.mass_model = mass_model_path

		if not config_revan == None:
			self.config_revan = config_revan

		if not config_mimrec == None:
			self.config_mimrec = config_mimrec

		if 'input_path' in inputs.keys() and 'output_path' in inputs.keys():
			if 'mass_model_path' in inputs.keys() and 'config_file_path' in inputs.keys():
				[input_path, 
		 	 	 output_path, 
		 	 	 mass_model_path, 
		 	 	 config_file_path] = define_paths([inputs['input_path'], inputs['output_path'], inputs['mass_model_path'], inputs['config_file_path']], 
		 									  	  [False, True, False, False])
				mass_model_path = mass_model_path[:-1]
				config_file_path = config_file_path[:-1]
				if not hasattr(self, 'mass_model'):
					self.mass_model = mass_model_path
				if not hasattr(self, 'config_revan'):
					self.config_revan = config_file_path
				if not hasattr(self, 'config_mimrec'):
					self.config_mimrec = config_file_path
			else:
				[input_path, output_path] = define_paths([inputs['input_path'], inputs['output_path']], [False, True])
			
			if not hasattr(self, 'input_cosima'):
				self.input_cosima = input_path
				self.input_revan = input_path
				self.input_mimrec = input_path
			if not hasattr(self, 'output_cosima'):
				self.output_cosima = output_path
				self.output_revan = output_path
				self.output_mimrec = output_path

	def cosima(self, file, output_path):
		"""
		Run cosima on a source file.

		Parameters
		----------
		file : str
			Source file name 
		"""

		if not os.path.isdir(output_path + 'output/'):
					os.makedirs(output_path + 'output/')
		os.system('cosima -z ' + file + ' > ' + output_path + 'output/output_' + file.split('.')[0] + '.txt')

	def mcosima(self, file, output_path, instances=None):
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
			os.system('mcosima -z ' + file + ' > ' + output_path + 'output/output_' + file.split('.')[0] + '.txt')
		else:
			os.system('mcosima -t ' + instances + ' -z ' + file + ' > ' + output_path + '/output_' + file.split('.')[0] + '.txt')

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
		os.system('mimrec -f ' + file + ' -g ' + mass_model + ' -c ' + config_mimrec + ' -n -x > ' + output_mimrec + 'output/output_' + file.split('.')[0] + '.txt')

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

	def run_cosima(self, input_path=None, output_path=None, zipped=True):
		"""
		Run cosima on all source files in a directory.

		Parameters
		----------
		input_path : str
			Path to directory housing source files
		output_path : str, optional
			Path to directory to store output sim or sim.gz files, if different from path defined when creating run_megalib object
		zipped : bool, optional
			Whether to zip output .sim files
		"""

		if input_path == None:
			input_path = self.input_cosima
		if output_path == None:
			output_path = self.output_cosima

		cosima_func = lambda file: self.cosima(file, output_path)

		if zipped:
			self.run_on_all(cosima_func, input_path, output_path, ['source'], 'sim.gz', 'Simulating')
		else:
			self.run_on_all(cosima_func, input_path, output_path, ['source'], 'sim', 'Simulating')

	def run_mcosima(self, input_path=None, output_path=None, instances=None, zipped=True):
		"""
		Run mcosima on all source files in a directory.

		Parameters
		----------
		input_path : str
			Path to directory housing source files
		output_path : str, optional
			Path to directory to store output sim or sim.gz files, if different from path defined when creating run_megalib object
		instances : int, optional
			Number of instances to run, default is number of CPU cores
		zipped : bool, optional
			Whether to zip output .sim files
		"""

		if input_path == None:
			input_path = self.input_cosima
		if output_path == None:
			output_path = self.output_cosima

		mcosima_func = lambda file: self.mcosima(file, output_path, instances=instances)

		if zipped:
			self.run_on_all(mcosima_func, input_path, output_path, ['source'], 'sim.gz', 'Simulating')
		else:
			self.run_on_all(mcosima_func, input_path, output_path, ['source'], 'sim', 'Simulating')

	def run_revan(self, input_path=None, output_path=None, config_file=None, zipped=True):
		"""
		Run revan on all sim and sim.gz files in a directory.

		Parameters
		----------
		input_path : str
			Path to directory housing sim and/or sim.gz files
		output_path : str, optional
			Path to directory to store output tra or tra.gz files, if different from path defined when creating run_megalib object
		zipped : bool, optional
			Whether to zip output tra files
		"""

		if input_path == None:
			input_path = self.input_revan
		if output_path == None:
			output_path = self.output_revan
		if config_file == None:
			config_file = self.config_revan

		revan_func = lambda file: self.revan(file, output_path, self.mass_model, config_file)

		if zipped:
			self.run_on_all(revan_func, input_path, output_path, ['sim','sim.gz'], 'tra.gz', 'Reconstructing')
		else:
			self.run_on_all(revan_func, input_path, output_path, ['sim','sim.gz'], 'tra', 'Reconstructing')

	def run_mimrec(self, input_path=None, output_path=None, config_file=None):
		"""
		Run mimrec on all tra and tra.gz files in a directory.

		Parameters
		----------
		input_path : str
			Path to directory housing tra and/or tra.gz files
		output_path : str, optional
			Path to directory to store output extracted tra files, if different from path defined when creating run_megalib object
		"""

		if input_path == None:
			input_path = self.input_mimrec
		if output_path == None:
			output_path = self.output_mimrec
		if config_file == None:
			config_file = self.config_mimrec

		mimrec_func = lambda file: self.mimrec(file, output_path, self.mass_model, config_file)

		self.run_on_all(mimrec_func, input_path, output_path, ['tra','tra.gz'], 'extracted.tra', 'Extracting events for')