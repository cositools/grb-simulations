import ROOT as root
import os
import random
import shutil
import pandas as pd
import h5py
import csv
import numpy as np
from scipy import integrate
from .config import suppress_output, read_yaml, define_paths
from .load_megalib import load_megalib
from .model import model

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

		self.config = inputs['configuration']

		if self.config == 'dc2':
			[self.source_path, 
		 	 self.background_path, 
		 	 self.output_path,
		 	 self.source_file_path, 
		 	 self.mass_model_path] = define_paths([inputs['paths']['input'], inputs['paths']['background'], inputs['paths']['output'], inputs['paths']['source_files'], inputs['paths']['mass_model']], 
		 									  	  [False, False, True, False, False])

			init_dc2_mass_model(inputs['mass_model_version'])

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

			if hasattr(self, 'background_length') and self.background_length >= 86400:
				raise RuntimeError('For now, background simulations must be shorter than a day to avoid potential error with make_hit_dict()')

		elif self.config == 'dc3':

			if 'background' in inputs['paths'] and 'saa_intervals' in inputs['paths'] and 'gbm_data' in inputs['paths']:
				[self.source_path, 
		 	 	 self.background_path, 
		 	 	 self.output_path,
		 	 	 self.source_file_path, 
		 	 	 self.mass_model_path,
		 	 	 self.saa_interval_path,
		 	 	 self.gbm_data_path] = define_paths([inputs['paths']['input'], inputs['paths']['background'], inputs['paths']['output'], inputs['paths']['source_files'], inputs['paths']['mass_model'], inputs['paths']['saa_intervals'], inputs['paths']['gbm_data']], 
		 												[False, False, True, False, False, False, False])

				___, self.saa_intervals = self.read_time_intervals(self.saa_interval_path)

			elif 'background' in inputs['paths'] and 'saa_intervals' in inputs['paths']:
				[self.source_path, 
		 	 	 self.background_path, 
		 	 	 self.output_path,
		 	 	 self.source_file_path, 
		 	 	 self.mass_model_path,
		 	 	 self.saa_interval_path] = define_paths([inputs['paths']['input'], inputs['paths']['background'], inputs['paths']['output'], inputs['paths']['source_files'], inputs['paths']['mass_model'], inputs['paths']['saa_intervals']], 
		 												[False, False, True, False, False, False])

				___, self.saa_intervals = self.read_time_intervals(self.saa_interval_path)

		 	elif 'background' in inputs['paths']:
				[self.source_path, 
		 	 	 self.background_path, 
		 	 	 self.output_path,
		 	 	 self.source_file_path, 
		 	 	 self.mass_model_path] = define_paths([inputs['paths']['input'], inputs['paths']['background'], inputs['paths']['output'], inputs['paths']['source_files'], inputs['paths']['mass_model']], 
		 									  		  [False, False, True, False, False])

			else:
		 		[self.source_path,  
		 		 self.output_path,
		 		 self.source_file_path, 
		 		 self.mass_model_path] = define_paths([inputs['paths']['input'], inputs['paths']['output'], inputs['paths']['source_files'], inputs['paths']['mass_model']], 
		 									  		  [False, True, False, False])

			self.event_list = self.read_event_list(self.source_file_path + 'event_list.txt')

			bgob1_pos = [3., 7.]
			bgob2_pos = [3., 7.]
			bgox1_pos = [15., -10.]
			bgox2_pos = [-10., -10.]
			bgoy1_pos = [15., 7.]
			bgoy2_pos = [-10., 7.]

			self.bgo_pos = [bgob1_pos, bgob2_pos, bgox1_pos, bgox2_pos, bgoy1_pos, bgoy2_pos]
			self.detector_keys = ['bgox1', 'bgox2', 'bgoy1', 'bgoy2', 'bgob1', 'bgob2']

		else:
			raise RuntimeError('Configuration must be dc2 or dc3')

		self.megalib = load_megalib(self.mass_model_path)

	def read_time_intervals(self, file_path):
		"""
		Read file with time intervals.

		Parameters
		----------
		file_path : str
			Path to file
		"""

		time_intervals = []

		with open(file_path, 'r') as csvfile:
			reader = csv.reader(csvfile)
			columns = next(reader)
			for row in reader:
				time_intervals.append(tuple([float(row[0]), float(row[1])]))

		return columns, time_intervals

	def init_dc2_mass_model(self, mass_model_version):
		"""
		Initialize mass model variables.

		Parameters
		----------
		mass_model_version : str
			Version of DC2 mass model ('8' or '12')
		"""

		if mass_model_version == 8:
			bgob_pos = [-19.2, 19.2, -17.5, 17.5, 14.6, 16.5]
			bgox1_pos = [17.5, 19.3, -17.5, 17.5, 16.9, 34.7]
			bgox2_pos = [-19.3, -17.5, -17.5, 17.5, 16.9, 34.7]
			bgoy1_pos = [-17.1, 17.1, 15.8, 17.6, 16.9, 34.7]
			bgoy2_pos = [-17.1, 17.1, -17.6, -15.8, 16.9, 34.7]
			ged1_z = [23.1, 24.5]
			ged2_z = [25.6, 27.0]
			ged3_z = [28.2, 29.6]
			ged4_z = [30.7, 32.1]

			self.bgo_pos = [bgob_pos, bgox1_pos, bgox2_pos, bgoy1_pos, bgoy2_pos]
			self.ged_pos = [ged1_z, ged2_z, ged3_z, ged4_z]
			self.detector_keys = ['ged1', 'ged2', 'ged3', 'ged4', 'bgox1', 'bgox2', 'bgoy1', 'bgoy2', 'bgob']

		elif mass_model_version == 12:
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

		else:
			raise RuntimeError('Mass model version ' + str(mass_model_version) + ' is not supported. Must be 8 or 12 for dc2')

	def read_event_list(self, event_list):
		"""
		Read event list file and save as a dictionary.

		Parameters
		----------
		event_list : dict
			Event list data
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

	def read_events_dc2(self, times, energies, hit):
		'''
		Sort hits into detectors based on positions.

		Parameters
		----------
		times : dict
			Times of hits in each detector
		energies : dict
			Energies of hits in each detector
		hit : ?
			New hit to add to dictionaries

		Returns
		----------
		times : dict
			Times of hits in each detector with new hit added
		energies : dict
			Energies of hits in each detector with new hit added
		'''

		if hit.GetDetector() == self.bgo_num and energy >= self.bgo_elim:
			if len(self.bgo_pos) == 5:
				# This may not be correct, events missing
				if self.bgo_pos[0][0] <= position[0] <= self.bgo_pos[0][1] and self.bgo_pos[0][2] <= position[1] <= self.bgo_pos[0][3] and self.bgo_pos[0][4] <= position[2] <= self.bgo_pos[0][5]:
					times['bgob'].append(float(time))
					energies['bgob'].append(float(energy))
				elif self.bgo_pos[1][0] <= position[0] <= self.bgo_pos[1][1] and self.bgo_pos[1][2] <= position[1] <= self.bgo_pos[1][3] and self.bgo_pos[1][4] <= position[2] <= self.bgo_pos[1][5]:
					times['bgox1'].append(float(time))
					energies['bgox1'].append(float(energy))
				elif self.bgo_pos[2][0] <= position[0] <= self.bgo_pos[2][1] and self.bgo_pos[2][2] <= position[1] <= self.bgo_pos[2][3] and self.bgo_pos[2][4] <= position[2] <= self.bgo_pos[2][5]:
					times['bgox2'].append(float(time))
					energies['bgox2'].append(float(energy))
				elif self.bgo_pos[3][0] <= position[0] <= self.bgo_pos[3][1] and self.bgo_pos[3][2] <= position[1] <= self.bgo_pos[3][3] and self.bgo_pos[3][4] <= position[2] <= self.bgo_pos[3][5]:
					times['bgoy1'].append(float(time))
					energies['bgoy1'].append(float(energy))
				elif self.bgo_pos[4][0] <= position[0] <= self.bgo_pos[4][1] and self.bgo_pos[4][2] <= position[1] <= self.bgo_pos[4][3] and self.bgo_pos[4][4] <= position[2] <= self.bgo_pos[4][5]:
					times['bgoy2'].append(float(time))
					energies['bgoy2'].append(float(energy))

		elif hit.GetDetector() == self.ged_num and hasattr(self, 'ged_pos'):
			if self.ged_pos[0][0] <= position[2] <= self.ged_pos[0][1]:
				times['ged1'].append(float(time))
				energies['ged1'].append(float(energy))
			elif self.ged_pos[1][0] <= position[2] <= self.ged_pos[1][1]:
				times['ged2'].append(float(time))
				energies['ged2'].append(float(energy))
			elif self.ged_pos[2][0] <= position[2] <= self.ged_pos[2][1]:
				times['ged3'].append(float(time))
				energies['ged3'].append(float(energy))
			elif self.ged_pos[3][0] <= position[2] <= self.ged_pos[3][1]:
				times['ged4'].append(float(time))
				energies['ged4'].append(float(energy))

		return times, energies

	def read_events_dc3(self, times, energies, event):
		'''
		Sort hits into detectors based on positions.

		Parameters
		----------
		times : dict
			Times of hits in each detector
		energies : dict
			Energies of hits in each detector
		event : ?
			New event to add to dictionaries

		Returns
		----------
		times : dict
			Times of hits in each detector with new event added
		energies : dict
			Energies of hits in each detector with new event added
		'''

		bottom_Zplus_1 = 0.
		bottom_Zplus_2 = 0.
		bottom_Zplus_3 = 0.
		bottom_Zplus_4 = 0.
		bottom_Zplus_5 = 0.
		bottom_Zminus_1 = 0.
		bottom_Zminus_2 = 0.
		bottom_Zminus_3 = 0.
		bottom_Zminus_4 = 0.
		bottom_Zminus_5 = 0.
		x1_1 = 0.
		x1_2 = 0.
		x1_3 = 0.
		x2_1 = 0.
		x2_2 = 0.
		x2_3 = 0.
		y1_1 = 0.
		y1_2 = 0.
		y1_3 = 0.
		y2_1 = 0.
		y2_2 = 0.
		y2_3 = 0.

		time = float(event.GetTime().GetAsSeconds())

		for i in range(event.GetNHTs()):
			hit = event.GetHTAt(i)

			if hit.GetDetectorType() == 8:
				x=hit.GetPosition().X()
				y=hit.GetPosition().Y()
				z=hit.GetPosition().Z()
		
				#bottom_Zplus
				if x<-11.2 and y>3 and z<7:
					bottom_Zplus_1+=(hit.GetEnergy())
				elif x>-11.2 and x<-2.5 and y>3 and z<7:
					bottom_Zplus_2+=(hit.GetEnergy())
				elif x>-2.5 and x<4.6 and y>3 and z<7:
					bottom_Zplus_3+=(hit.GetEnergy())
				elif x>4.6 and x<11.2 and y>3 and z<7:
					bottom_Zplus_4+=(hit.GetEnergy())
				elif x>4.6 and y>3 and z<7:
					bottom_Zplus_5+=(hit.GetEnergy())
               
				#bottom Zminus
				elif x<-11.2 and y<3 and z<7:
					bottom_Zminus_1+=(hit.GetEnergy())
				elif x>-11.2 and x<-2.5 and y<3 and z<7:
					bottom_Zminus_2+=(hit.GetEnergy())
				elif x>-2.5 and x<4.6 and y<3 and z<7:
					bottom_Zminus_3+=(hit.GetEnergy())
				elif x>4.6 and x<11.2 and y<3 and z<7:
 					bottom_Zminus_4+=(hit.GetEnergy())
				elif x>4.6 and y<3 and z<7:
					bottom_Zminus_5+=(hit.GetEnergy())
                
				#y1 pannel
				elif z>7 and y <-2.6 and y>-14 and x>15:
					y1_1+=(hit.GetEnergy())
				elif z>7 and y >-2.6 and y <9 and x>15:
					y1_2+=(hit.GetEnergy())
				elif z>7 and y >9 and y<20.6 and x>15:
					y1_3+=(hit.GetEnergy())
                
				#y2 pannel
				elif z>7 and y <-2.6 and y>-14 and x<-10:
					y2_1+=(hit.GetEnergy())
				elif z>7 and y >-2.6 and y <9 and x<-10:
					y2_2+=(hit.GetEnergy())
				elif z>7 and y >9 and y<20.6 and  x<-10:
					y2_3+=(hit.GetEnergy())

				#x1 panel
				elif z>-10 and x<-6 and y>15:
					x1_1+=(hit.GetEnergy())
				elif z>-10 and x >-6 and x <6 and y>15:
					x1_2+=(hit.GetEnergy())
				elif z>-10 and x >6 and y>15:
					x1_3+=(hit.GetEnergy())

				#x2 panel
				elif z>-10 and x<-6 and y<-10:
					x2_1+=(hit.GetEnergy())
				elif z>-10 and x >-6 and x <6 and y<-10:
					x2_2+=(hit.GetEnergy())
				elif z>-10 and x >6 and y<-10:
					x2_3+=(hit.GetEnergy())
				
				else:
					print(str(x)+" "+str(y)+" "+str(z))
					print("Error: Coordinate not found")

		if bottom_Zplus_1 >= 80.:
			times['bgob1'].append(float(time))
			energies['bgob1'].append(float(bottom_Zplus_1))
		if bottom_Zplus_2 >= 80.:
			times['bgob1'].append(float(time))
			energies['bgob1'].append(float(bottom_Zplus_2))
		if bottom_Zplus_3 >= 80.:
			times['bgob1'].append(float(time))
			energies['bgob1'].append(float(bottom_Zplus_3))
		if bottom_Zplus_4 >= 80.:
			times['bgob1'].append(float(time))
			energies['bgob1'].append(float(bottom_Zplus_4))
		if bottom_Zplus_5 >= 80.:
			times['bgob1'].append(float(time))
			energies['bgob1'].append(float(bottom_Zplus_5))
		if bottom_Zminus_1 >= 80.:
			times['bgob2'].append(float(time))
			energies['bgob2'].append(float(bottom_Zminus_1))
		if bottom_Zminus_2 >= 80.:
			times['bgob2'].append(float(time))
			energies['bgob2'].append(float(bottom_Zminus_2))
		if bottom_Zminus_3 >= 80.:
			times['bgob2'].append(float(time))
			energies['bgob2'].append(float(bottom_Zminus_3))
		if bottom_Zminus_4 >= 80.:
			times['bgob2'].append(float(time))
			energies['bgob2'].append(float(bottom_Zminus_4))
		if bottom_Zminus_5 >= 80.:
			times['bgob2'].append(float(time))
			energies['bgob2'].append(float(bottom_Zminus_5))
		if x1_1 >= 80.:
			times['bgox1'].append(float(time))
			energies['bgox1'].append(float(x1_1))
		if x1_2 >= 80.:
			times['bgox1'].append(float(time))
			energies['bgox1'].append(float(x1_2))
		if x1_3 >= 80.:
			times['bgox1'].append(float(time))
			energies['bgox1'].append(float(x1_3))
		if x2_1 >= 80.:
			times['bgox2'].append(float(time))
			energies['bgox2'].append(float(x2_1))
		if x2_2 >= 80.:
			times['bgox2'].append(float(time))
			energies['bgox2'].append(float(x2_2))
		if x2_3 >= 80.:
			times['bgox2'].append(float(time))
			energies['bgox2'].append(float(x2_3))
		if y1_1 >= 80.:
			times['bgoy1'].append(float(time))
			energies['bgoy1'].append(float(y1_1))
		if y1_2 >= 80.:
			times['bgoy1'].append(float(time))
			energies['bgoy1'].append(float(y1_2))
		if y1_3 >= 80.:
			times['bgoy1'].append(float(time))
			energies['bgoy1'].append(float(y1_3))
		if y2_1 >= 80.:
			times['bgoy2'].append(float(time))
			energies['bgoy2'].append(float(y2_1))
		if y2_2 >= 80.:
			times['bgoy2'].append(float(time))
			energies['bgoy2'].append(float(y2_2))
		if y2_3 >= 80.:
			times['bgoy2'].append(float(time))
			energies['bgoy2'].append(float(y2_3))

		return times, energies

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

		reach_end = False

		with suppress_output():
			while True: 
				event = reader.GetNextEvent()
				if not event or reach_end:
					break
				root.SetOwnership(event, True)
				time = float(event.GetTime().GetAsSeconds())

				if self.config == 'dc2':

					for i in range(event.GetNHTs()):
						hit = event.GetHTAt(i)
  
						if not end_time == None:
							if time > end_time:
								reach_end = True
								break
					times, energies = self.read_events_dc2(times, energies, hit)

				else:

					times, energies = self.read_events_dc3(times, energies, event)

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

	def read_background_sim(self, megalib, path, end_time):
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

	def read_background_csv(self, path, end_time=None):
		"""
		Fill hit directories for each background component.

		Parameters
		----------
		path : str
			Path to background file
		end_time : int or float
			Time at which to stop reading background .csv file

		Returns
		----------
		times : dict
			Times of hits in each detector
		energies : dict
			Energies of hits in each detector
		"""

		times = {'bgob1': [], 'bgob2': [], 'bgox1': [], 'bgox2': [], 'bgoy1': [], 'bgoy2': []}
		energies = {'bgob1': [], 'bgob2': [], 'bgox1': [], 'bgox2': [], 'bgoy1': [], 'bgoy2': []}

		print('Reading background file: ' + path.split('/')[-1])
		data = pd.read_csv(path, compression='gzip')
		
		for i in range(len(data['timestamp[s]'])):
			if not end_time == None:
				if data['timestamp[s]'][i] > end_time:
					break
			if data['bgo_bottom_1[keV]'][i] != 0.0:
				times['bgob1'].append(float(data['timestamp[s]'][i]))
				energies['bgob1'].append(float(data['bgo_bottom_1[keV]'][i]))
			elif data['bgo_bottom_2[keV]'][i] != 0.0:
				times['bgob2'].append(float(data['timestamp[s]'][i]))
				energies['bgob2'].append(float(data['bgo_bottom_2[keV]'][i]))
			elif data['bgo_x1[keV]'][i] != 0.0:
				times['bgox1'].append(float(data['timestamp[s]'][i]))
				energies['bgox1'].append(float(data['bgo_x1[keV]'][i]))
			elif data['bgo_x2[keV]'][i] != 0.0:
				times['bgox2'].append(float(data['timestamp[s]'][i]))
				energies['bgox2'].append(float(data['bgo_x2[keV]'][i]))
			elif data['bgo_y1[keV]'][i] != 0.0:
				times['bgoy1'].append(float(data['timestamp[s]'][i]))
				energies['bgoy1'].append(float(data['bgo_y1[keV]'][i]))
			elif data['bgo_y2[keV]'][i] != 0.0:
				times['bgoy2'].append(float(data['timestamp[s]'][i]))
				energies['bgoy2'].append(float(data['bgo_y2[keV]'][i]))

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
			if 'sim' in background_paths[i].split('.'):
				component_times, component_energies = self.read_background_sim(megalib, background_paths[i], end_time)
			elif 'csv' in background_paths[i].split('.'):
				component_times, component_energies = self.read_background_csv(background_paths[i], end_time)
			for key in component_times.keys():
				for j in range(len(component_times[key])):
					if component_times[key][j] > start_time:
						times[key].append(component_times[key][j])
						energies[key].append(component_energies[key][j])

		return times, energies

	def combine_single_source(self, source_times, source_energies, background_times, background_energies):
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

		data = np.array([times, energies])

		print('Writing to file: ' + file_path.split('/')[-1])
		with h5py.File(file_path, 'w') as f:
			dset = f.create_dataset('trigger_data', data=data, compression='gzip')
			dset.attrs['columns'] = ['time (s)', 'energy (keV)']

	def copy_readme(self):
		"""
		Copy README from source file directory to trigger input directory.
		"""

		if os.path.isfile(self.source_file_path + 'README.md'):
			shutil.copy(self.source_file_path + 'README.md', self.output_path + 'README.md')
		else:
			print('No README in .source file directory.')

	def write_readme_dc2(self, source_name, start_time):
		"""
		Write README for event directory.

		Parameters
		----------
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

	def write_event_readme_dc3(self, source_name):
		"""
		Write README for event directory.

		Parameters
		----------
		source_name : str
			Name of source
		"""

		source_index = np.where(np.array(self.event_list['Event name'])==source_name)[0][0]

		with open(self.output_path + source_name + '/' + 'event_list.txt', 'w') as f:
			f.write('zenith angle (deg): ' + self.event_list[' Azimuth (degrees)'] + '\n')
			f.write('azimuth angle (deg): ' + self.event_list[' Zenith (degrees)'] + '\n')
			f.write('energy flux (erg/cm^2/s): ' + self.event_list[' Energy flux (erg/cm^2/s)'] + '\n')
			f.write('photon flux (ph/cm^2/s): ' + self.event_list[' Photon flux (ph/cm^2/s)'] + '\n')
			f.write('event start time (s): ' + str(self.event_list[' Start time (s)']) + '\n')
			f.write('event duration (s): ' + str(self.event_list[' Duration (s)']) + '\n')

	def write_readme_dc3(self, directory_path, event_list):
		"""
		Write README for directory of multiple events combined with background.

		Parameters
		----------
		directory_path : str
			Path to directory containing trigger files
		event_list : list
			Event list for events in directory
		"""

		with open(directory_path + 'event_list.txt', 'w') as f:

			if not hasattr(self, 'gbm_data_path'):

				f.write('Event name, Start time (s), Duration (s), Photon flux (ph/cm^2/s), Energy flux (erg/cm^2/s), Zenith (degrees), Azimuth (degrees)\n')
				for i in range(len(event_list['Event name'])):
					f.write(event_list['Event name'][i] + ', ' + str(event_list['Start time (s)'][i]) + ', ' + str(event_list['Duration (s)'][i]) + ', ' + str(event_list['Photon flux (ph/cm^2/s)'][i]) + ', ' + 
							str(event_list['Energy flux (erg/cm^2/s)'][i]) + ', ' + str(event_list['Zenith (degrees)'][i]) + ', ' + str(event_list['Azimuth (degrees)'][i]) + '\n')

			else:

				f.write('Event name, Start time (s), Duration (s), Photon flux (ph/cm^2/s) for 80-2000 keV, Photon flux (ph/cm^2/s) for 10-10000 keV, Energy flux (erg/cm^2/s) for 10-10000 keV, Zenith (degrees), Azimuth (degrees), Spectral model, Spectral parameters, GBM T90 (s)\n')
				for i in range(len(event_list['Event name'])):

					name = event_list['Event name'][i]
					event_file = read_yaml(self.gbm_data_path + name + '/' + name + '.yaml')
					spectral_model = spectrum_file['flnc_best_fitting_model']

					if spectral_model == 'flnc_comp':
						spectral_model = 'Comptonized'
						parameters = [event_file['flnc_comp_index'], event_file['flnc_comp_epeak']]
						parameters_file = parameters
					elif spectral_model == 'flnc_band':
						spectral_model = 'Band'
						parameters = [event_file['flnc_band_alpha'], this_file['event_band_beta'], this_file['flnc_band_epeak'] / (this_file['flnc_band_alpha'] + 2)]
						parameters_file = [event_file['flnc_band_alpha'], this_file['event_band_beta'], this_file['event_band_epeak']]
					else:
						raise RuntimeError(spectral_model + ' spectral model not supported')

					this_model = model(spectral_model)
					this_model.set_model(parameters)
					bgo_flux = float(event_list['Photon flux (ph/cm^2/s)'][i]) / integrate.quad(this_model.function, 10, 10000)[0] * integrate.quad(this_model.function, 80, 2000)[0]

					f.write(event_list['Event name'][i] + ', ' + str(event_list['Start time (s)'][i]) + ', ' + str(event_list['Duration (s)'][i]) + ', ' + f'{bgo_flux:.4f}' + ', ' + str(event_list['Photon flux (ph/cm^2/s)'][i]) + ', ' + 
							str(event_list['Energy flux (erg/cm^2/s)'][i]) + ', ' + str(event_list['Zenith (degrees)'][i]) + ', ' + str(event_list['Azimuth (degrees)'][i]) + ', ' + spectral_model + ', ' + parameters_file + ', ' + str(event_file['t90']) + '\n')

	def combine_multiple_events(self, source_times, source_energies, background_times, background_energies):
		"""
		Combine source & background events and sort by time.

		Parameters
		----------
		source_times : dict
			Times of source events
		source_energies : dict
			Energies of source events
		background_times : dict
			Times of source events
		background_energies : dict
			Energies of source events

		Returns
		----------
		times_sorted : dict
			Times of hits in each detector
		energies_sorted : dict
			Energies of hits in each detector
		"""

		for key in background_times.keys():
			for i in range(len(background_times[key])):
				source_times[key].append(background_times[key][i])
		for key in background_energies.keys():
			for i in range(len(background_energies[key])):
				source_energies[key].append(background_energies[key][i])
		times_sorted = {}
		energies_sorted = {}
		for item in self.detector_keys:
			times_sorted[item] = []
			energies_sorted[item] = []
		for key in background_times.keys():
			times_sorted[key], energies_sorted[key] = (list(x) for x in zip(*sorted(zip(source_times[key], source_energies[key]))))

		return times_sorted, energies_sorted

	def create_event_files(self):
		"""
		Create trigger algorithm input files for all events in input directory.
		"""

		self.copy_readme()

		if self.config == 'dc3':

			if hasattr(self, 'background_path'):
				n_background_files = 0
				for background_file in os.listdir(self.background_path):
					if n_background_files == 0:
						background_times, background_energies = self.read_background_csv(self.background_path + background_file)
						break
					else:
						component_times, component_energies = self.read_background_csv(self.background_path + background_file)
						for key in background_times.keys():
							background_times[key].extend(component_times[key])
							background_energies[key].extend(component_energies[key])
					n_background_files += 1
				for key in background_times.keys():
					background_times[key], background_energies[key] = (list(x) for x in zip(*sorted(zip(background_times[key], background_energies[key]))))

				directory_number = 0
				previous_time = -1.
				for i in range(len(self.event_list['Event name'])):

					this_time = float(self.event_list[' Start time (s)'][i])

					if previous_time > this_time:
						directory_path = self.output_path + 'batch_' + str(directory_number) + '/'
						os.mkdir(directory_path)
						self.write_readme_dc3(directory_path, this_event_list)
						times_sorted, energies_sorted = self.combine_multiple_events(batch_times, batch_energies, background_times, background_energies)
						for key in times_sorted.keys():
							self.write_hits(directory_path + key + '.hdf5', times_sorted[key], energies_sorted[key])
						directory_number += 1

					filename = None
					for file in os.listdir(self.source_path):
						if file.startswith(self.event_list['Event name'][i]):
							filename = file
					if filename is None:
						print('.sim file for ' + self.event_list['Event name'][i] + ' not found')
						if previous_time > this_time or i == 0:
							this_event_list= {'Event name': [], 'Start time (s)': [], 'Duration (s)': [], 'Photon flux (ph/cm^2/s)': [], 'Energy flux (erg/cm^2/s)': [], 'Zenith (degrees)': [], 'Azimuth (degrees)': []}
							batch_times = {}
							batch_energies = {}
							for key in self.detector_keys:
								batch_times[key] = []
								batch_energies[key] = []
						if i == len(self.event_list['Event name']) - 1:
							directory_path = self.output_path + 'batch_' + str(directory_number) + '/'
							os.mkdir(directory_path)
							self.write_readme_dc3(directory_path, this_event_list)
							times_sorted, energies_sorted = self.combine_multiple_events(batch_times, batch_energies, background_times, background_energies)
							for key in times_sorted.keys():
								self.write_hits(directory_path + key + '.hdf5', times_sorted[key], energies_sorted[key])
						continue

					event_list_values = []
					for key in self.event_list:
						event_list_values.append(self.event_list[key][i])
    					
					self.megalib.open_file(self.source_path + filename)

					if previous_time > this_time or i == 0:

						this_event_list= {'Event name': [], 'Start time (s)': [], 'Duration (s)': [], 'Photon flux (ph/cm^2/s)': [], 'Energy flux (erg/cm^2/s)': [], 'Zenith (degrees)': [], 'Azimuth (degrees)': []}
						saa = False
						if hasattr(self, 'saa_intervals'):
							end_time = this_time + float(self.event_list[' Duration (s)'][i])
							for interval in self.saa_intervals:
								if (this_time >= interval[0] and this_time <= interval[1]) or (end_time >= interval[0] and end_time <= interval[1]):
									print(filename + ' occurs during SAA passage and will not be included')
									batch_times = {}
									batch_energies = {}
									for key in self.detector_keys:
										batch_times[key] = []
										batch_energies[key] = []
									saa = True
									break
						if not saa:
							for key, value in zip(this_event_list, event_list_values):
								this_event_list[key].append(value)
							print('Reading sim file: ' + filename)
							batch_times, batch_energies = self.make_hit_dict(self.megalib.reader)
				
					else:

						saa = False
						if hasattr(self, 'saa_intervals'):
							end_time = this_time + float(self.event_list[' Duration (s)'][i])
							for interval in self.saa_intervals:
								if (this_time >= interval[0] and this_time <= interval[1]) or (end_time >= interval[0] and end_time <= interval[1]):
									print(filename + ' occurs during SAA passage and will not be included')
									saa = True
									break
						if not saa:
							for key, value in zip(this_event_list, event_list_values):
								this_event_list[key].append(value)
							print('Reading sim file: ' + filename)
							these_times, these_energies = self.make_hit_dict(self.megalib.reader)
							for key in these_times.keys():
								for item in these_times[key]:
									batch_times[key].append(item)
							for key in these_energies.keys():
								for item in these_energies[key]:
									batch_energies[key].append(item)

						if i == len(self.event_list['Event name']) - 1:
							directory_path = self.output_path + 'batch_' + str(directory_number) + '/'
							os.mkdir(directory_path)
							self.write_readme_dc3(directory_path, this_event_list)
							times_sorted, energies_sorted = self.combine_multiple_events(batch_times, batch_energies, background_times, background_energies)
							for key in times_sorted.keys():
								self.write_hits(directory_path + key + '.hdf5', times_sorted[key], energies_sorted[key])
					
					previous_time = this_time

			else:

				for file in os.listdir(self.source_path):
					filename = os.fsdecode(file)
					source_name = file.split('.')[0]

					if not os.path.isdir(self.output_path + source_name) and os.path.isfile(self.source_path + file) and (file[-3:] == 'sim' or file[-6:] == 'sim.gz'):
						self.megalib.open_file(self.source_path + filename)
						os.mkdir(self.output_path + source_name + '/')
						print('Reading sim file: ' + filename)
						source_times, source_energies = self.make_hit_dict(self.megalib.reader)

						for key in times.keys():
							self.write_hits(self.output_path + source_name + '/' + key + '.hdf5', times[key], energies[key])

						self.write_event_readme_dc3(source_name)

		else:
			for file in os.listdir(self.source_path):
				filename = os.fsdecode(file)
				source_name = file.split('.')[0]

				if not os.path.isdir(self.output_path + source_name) and os.path.isfile(self.source_path + file) and (file[-3:] == 'sim' or file[-6:] == 'sim.gz'):
					self.megalib.open_file(self.source_path + filename)
					os.mkdir(self.output_path + source_name + '/')
					print('Reading sim file: ' + filename)
					source_times, source_energies = self.make_hit_dict(self.megalib.reader)
    
					if self.background_type == 'random':
						background_paths, start, end = self.choose_background()
						background_times, background_energies = self.select_background(self.megalib, background_paths, start, end)
					elif self.background_type == 'file':
						self.megalib.open_file(self.background_path)
						print('Reading background file: ' + self.background_path.split('/')[-1])
						background_times, background_energies = self.make_hit_dict(self.megalib.reader)
					
					times, energies, start_time = self.combine_single_source(source_times, source_energies, background_times, background_energies)

					for key in times.keys():
						self.write_hits(self.output_path + source_name + '/' + key + '.hdf5', times[key], energies[key])

					self.write_readme_dc2(source_name, start_time)