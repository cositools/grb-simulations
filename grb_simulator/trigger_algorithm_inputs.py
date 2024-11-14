import ROOT as root
import os
import random
import shutil
from .config import suppress_output, read_yaml, define_paths
from .load_megalib import load_megalib

class trigger_algorithm_inputs():

	def __init__(self, input_file):
		"""
		Create trigger algorithm input files.

		Parameters
		----------
		input_file : str
			Path to input .yaml file
		"""

		inputs = read_yaml(input_file)

		[self.source_path, 
		 self.background_path, 
		 self.output_path,
		 self.source_file_path, 
		 self.mass_model_path] = define_paths([inputs['paths']['input'], inputs['paths']['background'], inputs['paths']['output'], inputs['paths']['source_files'], inputs['paths']['mass_model']], 
		 									  [False, False, True, False, False])
		if inputs['background']['type'] == 'file':
			self.background_path = self.background_path[:-1]
		self.mass_model_path = self.mass_model_path[:-1]

		if inputs['mass_model_version'] == 8 or inputs['mass_model_version'] == 12:
			self.mass_model_version = inputs['mass_model_version']
		else:
			raise RuntimeError('Mass model version ' + str(self.mass_model_version) + ' is not supported. Must be 8 or 12')

		if 'time' in inputs['background']:
			self.background_time = inputs['background']['time']
		if 'number' in inputs['background']:
			self.background_num = inputs['background']['number']
		if 'file_length' in inputs['background']:
			self.background_length = inputs['background']['file_length']
		if 'components' in inputs['background']:
			self.background_components = inputs['background']['components']
		if 'file_type' in inputs['background']:
			self.background_file_type = inputs['background']['file_type']
		if 'type' in inputs['background']:
			self.background_type = inputs['background']['type']

		if self.background_length >= 86400:
			raise RuntimeError('For now, background simulations must be shorter than a day to avoid potential error with make_hit_dict()')

		self.init_variables()

		self.megalib = load_megalib(self.mass_model_path)

	def init_variables(self):
		"""
		Initialize variables.
		"""

		self.bgo_num = 4
		self.ged_num = 3
		self.bgo_elim = 80

		if self.mass_model_version == 8:
			bgob_pos = [-19.2, 19.2, -17.5, 17.5, 14.6, 16.5]
			bgox1_pos = [17.5, 19.3, -17.5, 17.5, 16.9, 34.7]
			bgox2_pos = [-19.3, -17.5, -17.5, 17.5, 16.9, 34.7]
			bgoy1_pos = [-17.1, 17.1, 15.8, 17.6, 16.9, 34.7]
			bgoy2_pos = [-17.1, 17.1, -17.6, -15.8, 16.9, 34.7]
			ged1_z = [23.1, 24.5]
			ged2_z = [25.6, 27.0]
			ged3_z = [28.2, 29.6]
			ged4_z = [30.7, 32.1]

		elif self.mass_model_version == 12:
			bgob_pos = [-21.1, 21.1, -18.3, 18.3, 13.6, 16.]
			bgox1_pos = [19., 21.2, -18.3, 18.3, 16.3, 35.3]
			bgox2_pos = [-21.2, -19., -18.3, 18.3, 16.3, 35.3]
			bgoy1_pos = [-18.6, 18.6, 16.2, 18.4, 16.3, 35.3]
			bgoy2_pos = [-18.6, 18.6, -18.4, -16.2, 16.3, 35.3]
			ged1_z = [22.5, 23.9]
			ged2_z = [25.1, 26.5]
			ged3_z = [27.6, 29.0]
			ged4_z = [30.2, 31.6]

		self.bgo_pos = [bgob_pos, bgox1_pos, bgox2_pos, bgoy1_pos, bgoy2_pos]
		self.ged_pos = [ged1_z, ged2_z, ged3_z, ged4_z]

		self.detector_keys = ['ged1', 'ged2', 'ged3', 'ged4', 'bgox1', 'bgox2', 'bgoy1', 'bgoy2', 'bgob']

	def make_hit_dict(self, reader, end_time=None):
		"""
		Create and fill directories for each hit.

		Parameters
		----------
		reader : cppyy.gbl.MFileEventsSim
			MEGAlib reader
		end_time : int or float, optional
			Time at which to stop reading .sim file

		Returns
		----------
		times : dict
			Times of hits in each detector
		energies : dict
			Energies of hits in each detector
		"""
  
		times = {}
		energies = {}

		for item in self.detector_keys:
			times[item] = []
			energies[item] = []

		hrs_init = 19
		reach_end = False

		with suppress_output():
			while True: 
				event = reader.GetNextEvent()
				if not event or reach_end:
					break
				root.SetOwnership(event, True)

				for i in range(event.GetNHTs()):
					hit = event.GetHTAt(i)
					position = hit.GetPosition()
					energy = hit.GetEnergy()
					hour = int(event.GetTime().GetHours() - hrs_init)
					hour = int(hour + 24) if hour < 0 else hour
					nanosecond = str(event.GetTime().GetNanoSeconds()).zfill(9)
					seconds = str((hour * 3600) + (int(event.GetTime().GetMinutes()) * 60) + int(event.GetTime().GetSeconds()))
					time = float(seconds + '.' + nanosecond)
  
					if not end_time == None:
						if time > end_time:
							reach_end = True
							break
      
					if hit.GetDetector() == self.bgo_num and energy >= self.bgo_elim:
						if self.bgo_pos[0][0] <= position[0] <= self.bgo_pos[0][1] and self.bgo_pos[0][2] <= position[1] <= self.bgo_pos[0][3] and self.bgo_pos[0][4] <= position[2] <= self.bgo_pos[0][5]:
							times['bgob'].append(time)
							energies['bgob'].append(energy)
						elif self.bgo_pos[1][0] <= position[0] <= self.bgo_pos[1][1] and self.bgo_pos[1][2] <= position[1] <= self.bgo_pos[1][3] and self.bgo_pos[1][4] <= position[2] <= self.bgo_pos[1][5]:
							times['bgox1'].append(time)
							energies['bgox1'].append(energy)
						elif self.bgo_pos[2][0] <= position[0] <= self.bgo_pos[2][1] and self.bgo_pos[2][2] <= position[1] <= self.bgo_pos[2][3] and self.bgo_pos[2][4] <= position[2] <= self.bgo_pos[2][5]:
							times['bgox2'].append(time)
							energies['bgox2'].append(energy)
						elif self.bgo_pos[3][0] <= position[0] <= self.bgo_pos[3][1] and self.bgo_pos[3][2] <= position[1] <= self.bgo_pos[3][3] and self.bgo_pos[3][4] <= position[2] <= self.bgo_pos[3][5]:
							times['bgoy1'].append(time)
							energies['bgoy1'].append(energy)
						elif self.bgo_pos[4][0] <= position[0] <= self.bgo_pos[4][1] and self.bgo_pos[4][2] <= position[1] <= self.bgo_pos[4][3] and self.bgo_pos[4][4] <= position[2] <= self.bgo_pos[4][5]:
							times['bgoy2'].append(time)
							energies['bgoy2'].append(energy)
					elif hit.GetDetector() == self.ged_num:
						if self.ged_pos[0][0] <= position[2] <= self.ged_pos[0][1]:
							times['ged1'].append(time)
							energies['ged1'].append(energy)
						elif self.ged_pos[1][0] <= position[2] <= self.ged_pos[1][1]:
							times['ged2'].append(time)
							energies['ged2'].append(energy)
						elif self.ged_pos[2][0] <= position[2] <= self.ged_pos[2][1]:
							times['ged3'].append(time)
							energies['ged3'].append(energy)
						elif self.ged_pos[3][0] <= position[2] <= self.ged_pos[3][1]:
							times['ged4'].append(time)
							energies['ged4'].append(energy)
        
		return times, energies

	def choose_background(self):
		"""
		Randomly choose background files from specified directory.

		Returns
		----------
		background_paths : list of str
			Paths to background files
		background_start : int
			Time in background files to start saving hits
		background_end : int
			Time in background files to stop saving hits
		"""

		file_num = str(random.randint(1, self.background_num))

		if self.background_file_type == 'sequential':
			start_time = (int(file_num) - 1) * self.background_length
			background_start = random.randint(start_time, start_time + self.background_length - self.background_time)
			background_end = background_start + self.background_time
		elif self.background_file_type == 'simultaneous':
			background_start = random.randint(0, self.background_length - self.background_time)
			background_end = background_start + self.background_time
		else:
			raise RuntimeError("background_file_type in yaml file must be 'sequantial' or 'simultaneous'")

		background_paths = []

		for file in os.listdir(self.background_path):
			for item in self.background_components:
				if file.startswith(item + '.') and 'inc' + file_num + '.' in file and (file.split('.')[-1] == 'sim' or (file.split('.')[-1] == 'gz' and file.split('.')[-2] == 'sim')):
					background_paths.append(self.background_path + os.fsdecode(file))

		return background_paths, background_start, background_end

	def read_background(self, megalib, path, end_time):
		"""
		Fill hit directories for each background component.

		Parameters
		----------
		megalib : grb_simulator.load_megalib.load_megalib
			Load MEGAlib object
		path : str
			Path to background file
		end_time : int or float
			Time at which to stop reading background .sim file

		Returns
		----------
		times : dict
			Times of hits in each detector
		energies : dict
			Energies of hits in each detector
		"""

		megalib.open_file(path)
		print('Reading background file: ' + path.split('/')[-1])
		times, energies = self.make_hit_dict(megalib.reader, end_time)

		return times, energies

	def select_background(self, megalib, background_paths, start_time, end_time):
		"""
		Read in background events.

		Parameters
		----------
		megalib : grb_simulator.load_megalib.load_megalib
			Load MEGAlib object
		background_paths : list of str
			Paths to background files
		start_time : int or float
			Time at which to start saving background hits
		end_time : int or float
			Time at which to stop saving background hits

		Returns
		----------
		times : dict
			Times of hits in each detector
		energies : dict
			Energies of hits in each detector
		"""

		times = {}
		energies = {}

		for item in self.detector_keys:
			times[item] = []
			energies[item] = []

		for i in range(len(background_paths)):
			component_times, component_energies = self.read_background(megalib, background_paths[i], end_time)
			for key in component_times.keys():
				for j in range(len(component_times[key])):
					if component_times[key][j] > start_time:
						times[key].append(component_times[key][j])
						energies[key].append(component_energies[key][j])

		return times, energies

	def combine(self, source_times, source_energies, background_times, background_energies):
		"""
		Combine source and background hits by placing source in the middle of background.

		Parameters
		----------
		source_times : dict
			Times of hits from source in each detector
		source_energies : dict
			Energies of hits from source in each detector
		background_times : dict
			Times of hits from background in each detector
		background_energies : dict
			Energies of hits from background in each detector

		Returns
		----------
		times_sorted : dict
			Times of hits in each detector
		energies_sorted : dict
			Energies of hits in each detector
		source_start_time : float
			Time at which event begins
		"""

		min_times = []
		max_times = []
		times = {}
		energies = {}
		times_sorted = {}
		energies_sorted = {}

		for item in self.detector_keys:
			times[item] = []
			energies[item] = []
			times_sorted[item] = []
			energies_sorted[item] = []

		for key in background_times.keys():
			try:
				min_times.append(min(background_times[key]))
				max_times.append(max(background_times[key]))
			except:
				print(key, background_times[key])

		for key in source_times.keys():
			for i in range(len(source_times[key])):
				source_times[key][i] += (min(min_times) + max(max_times)) / 2

		source_start_time = min(source_times)

		for key in times.keys():
			for i in range(len(source_times[key])):
				times[key].append(source_times[key][i])
				energies[key].append(source_energies[key][i])
			for i in range(len(background_times[key])):
				times[key].append(background_times[key][i])
				energies[key].append(background_energies[key][i])
			times_sorted[key], energies_sorted[key] = (list(x) for x in zip(*sorted(zip(times[key], energies[key]))))

		return times_sorted, energies_sorted, source_start_time

	def write_hits(self, file_path, times, energies):
		"""
		Write times and energies to trigger algorithm input file.

		Parameters
		----------
		file_path : str
			Path to trigger algorithm input file
		times : dict
			Times of hits from source in each detector
		energies : dict
			Energies of hits from source in each detector
		"""

		print('Writing to file: ' + file_path.split('/')[-1])
		with open(file_path, 'w') as f:
			f.write('time (s)           energy (keV)\n')
			for i in range(len(times)):
				f.write("{:.9f}".format(times[i]) + '        ' + str(energies[i]) + '\n')

	def copy_readme(self):
		"""
		Copy README from source file directory to trigger input directory.
		"""

		if os.path.isfile(self.source_file_path + 'README.md'):
			shutil.copy(self.source_file_path + 'README.md', self.output_path + 'README.md')
		else:
			print('No README in .source file directory.')

	def write_readme(self, source_name, start_time):
		"""
		Write README for event directory.

		source_name : str
			Name of source
		start_time : float
			Time at which event begins
		"""

		filename = self.source_file_path + source_name + '.source'
		with open(filename, 'r') as f:
			lines = f.readlines()
		for line in lines:
			line_list = line.split()
			if '.Beam' in line:
				z = line_list[2]
				a = line_list[3]
			elif '.Spectrum' in line:
				spectrum = line_list[1:]
			elif 'Average photon flux' in line:
				e_flux = line_list[8:]
			elif '.Flux' in line:
				ph_flux = line_list[1] + ' ph/cm2/s'

		with open(self.output_path + source_name + '/' + 'README.md', 'w') as f:
			f.write('zenith angle: ' + z + '\n')
			f.write('azimuth angle: ' + a + '\n')
			f.write('spectrum: ' + spectrum[0] + '\n')
			f.write('spectral parameters or file: ' + str(spectrum[1:]) + '\n')
			f.write('energy flux: ' + str(e_flux[0]) + ' ' + str(e_flux[1]) + '\n')
			f.write('photon flux: ' + ph_flux + '\n')
			f.write('event start time: ' + str(start_time) + '\n')

	def create_event_files(self):
		"""
		Create trigger algorithm input files for all events in input directory.
		"""

		self.copy_readme()

		for file in os.listdir(self.source_path):
			filename = os.fsdecode(file)
			source_name = file.split('.')[0]
			if not os.path.isdir(self.output_path + source_name) and os.path.isfile(self.source_path + file) and (file[-3:] == 'sim' or file[-6:] == 'sim.gz'):
				self.megalib.open_file(self.source_path + filename)
				os.mkdir(self.output_path + source_name + '/')
				print('Reading source file: ' + filename)
				source_times, source_energies = self.make_hit_dict(self.megalib.reader)
    
				if self.background_type == 'random':
					background_paths, start, end = self.choose_background()
					background_times, background_energies = self.select_background(self.megalib, background_paths, start, end)
				else:
					self.megalib.open_file(self.background_path)
					print('Reading background file: ' + self.background_path.split('/')[-1])
					background_times, background_energies = self.make_hit_dict(self.megalib.reader)

				times, energies, start_time = self.combine(source_times, source_energies, background_times, background_energies)

				for key in times.keys():
					self.write_hits(self.output_path + source_name + '/' + key + '.txt', times[key], energies[key])

				self.write_readme(source_name, start_time)