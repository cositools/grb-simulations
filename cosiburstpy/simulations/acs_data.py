import shutil
import pandas as pd
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import logging
from cosiburstpy.utility.utility import read_hdf5, write_hdf5, SuppressOutput

logger = logging.getLogger(__name__)

class ACSData():

	def __init__(self, data):
		'''
		ACS data handling.

		Parameters
		----------
		data : dict
			ACS data where keys are ACS panel names and entries are lists of tuples of time and energy
		'''

		self.b1 = data['b1']
		self.b2 = data['b2']

		self.x1 = data['x1']
		self.x2 = data['x2']

		self.y1 = data['y1']
		self.y2 = data['y2']

	@classmethod
	def from_file(cls, file, mass_model=None):
		'''
		Read in ACS data file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to ACS data file
		mass_model : pathlib.PosixPath, optional
			Path to analysis mass model if reading .sim or .sim.gz file

		Returns
		-------
		acs_data : cosiburstpy.acs_data.ACSData
			ACS data
		'''

		if file.suffix == '.hdf5':

			acs_data = cls.from_hdf5_file(file)

		elif file.suffix == '.sim' or ''.join(file.suffixes) == '.sim.gz':

			acs_data = cls.from_sim_file(file, mass_model)

		elif ''.join(file.suffixes) == '.csv.gz':

			acs_data = cls.from_csv_file(file)

		else:

			raise RuntimeError(f'{file.suffix} files are not supported for ACS data.')

		return acs_data

	@classmethod
	def from_hdf5_file(cls, file):
		'''
		Extract ACS hits from .hdf5 file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to ACS data .hdf5 file
	
		Returns
		-------
		acs_data : cosiburstpy.acs_data.ACSData
			ACS data
		'''

		logger.info(f'Reading {file}')

		data, file_attributes, dataset_attributes = read_hdf5(file)
		data = {key: [tuple(value * unit for value, unit in zip(hit, (u.s, u.keV))) for hit in hits] for key, hits in data.items()}

		acs_data = cls(data)

		return acs_data

	@classmethod
	def from_sim_file(cls, file, mass_model):
		'''
		Extract ACS hits from .sim or .sim.gz file. Modified from Savitri's code to create ACS data .csv files with the Data Challenge 3 mass model.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to ACS data .sim or .sim.gz file
		mass_model : pathlib.PosixPath
			Path to analysis mass model
	
		Returns
		-------
		acs_data : cosiburstpy.acs_data.ACSData
			ACS data
		'''

		logger.info(f'Reading {file}')

		from cosiburstpy.megalib.read_sim_file import read_sim_file

		acs_data = read_sim_file(file, mass_model)

		return acs_data

	@classmethod
	def from_csv_file(cls, file):
		'''
		Extract ACS hits from .csv.gz file in Savitri's format.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to ACS data .csv.gz file
	
		Returns
		-------
		acs_data : cosiburstpy.acs_data.ACSData
			ACS data
		'''

		logger.info(f'Reading {file}')

		times = {'b1': [], 'b2': [], 'x1': [], 'x2': [], 'y1': [], 'y2': []}
		energies = {'b1': [], 'b2': [], 'x1': [], 'x2': [], 'y1': [], 'y2': []}

		data = pd.read_csv(file, compression='gzip')
		
		for i in range(len(data['timestamp[s]'])):

			if data['bgo_bottom_1[keV]'][i] != 0.0:

				times['b1'].append(float(data['timestamp[s]'][i]) * u.s)
				energies['b1'].append(float(data['bgo_bottom_1[keV]'][i]) * u.keV)

			elif data['bgo_bottom_2[keV]'][i] != 0.0:

				times['b2'].append(float(data['timestamp[s]'][i]) * u.s)
				energies['b2'].append(float(data['bgo_bottom_2[keV]'][i]) * u.keV)

			elif data['bgo_x1[keV]'][i] != 0.0:

				times['x1'].append(float(data['timestamp[s]'][i]) * u.s)
				energies['x1'].append(float(data['bgo_x1[keV]'][i]) * u.keV)

			elif data['bgo_x2[keV]'][i] != 0.0:

				times['x2'].append(float(data['timestamp[s]'][i]) * u.s)
				energies['x2'].append(float(data['bgo_x2[keV]'][i]) * u.keV)

			elif data['bgo_y1[keV]'][i] != 0.0:

				times['y1'].append(float(data['timestamp[s]'][i]) * u.s)
				energies['y1'].append(float(data['bgo_y1[keV]'][i]) * u.keV)

			elif data['bgo_y2[keV]'][i] != 0.0:

				times['y2'].append(float(data['timestamp[s]'][i]) * u.s)
				energies['y2'].append(float(data['bgo_y2[keV]'][i]) * u.keV)

		acs_data = cls({key: list(zip(times[key], energies[key])) for key in times})

		return acs_data

	@classmethod
	def combine(cls, files, mass_model=None):
		'''
		Combine ACS data files.

		Parameters
		----------
		files : list of pathlib.PosixPath
			Paths to ACS data files
		mass_model : pathlib.PosixPath, optional
			Path to analysis mass model if reading .sim file

		Returns
		-------
		acs_data : cosiburstpy.acs_data.ACSData
			Combined ACS data
		'''

		logger.info('Combining ACS data files')

		data = {'b1': [], 'b2': [], 'x1': [], 'x2': [], 'y1': [], 'y2': []}

		for file in files:

			component_data = cls.from_file(file, mass_model)

			data['b1'].extend(component_data.b1)
			data['b2'].extend(component_data.b2)
			data['x1'].extend(component_data.x1)
			data['x2'].extend(component_data.x2)
			data['y1'].extend(component_data.y1)
			data['y2'].extend(component_data.y2)

		for panel in data:
			data[panel].sort(key=lambda x: x[0].to(u.s).value)

		acs_data = cls(data)

		return acs_data

	@classmethod
	def convert(cls, input_file, output_file, mass_model=None):
		'''
		Convert ACS data file to .hdf5.

		Parameters
		----------
		input_file : pathlib.PosixPath
			Path to ACS data .sim, .sim.gz, or .csv.gz file
		output_file : pathlib.PosixPath
			Path to ACS data .hdf5 file
		mass_model : pathlib.PosixPath, optional
			Path to analysis mass model if reading .sim or .sim.gz file
		'''

		logger.info(f'Converting {input_file} to .hdf5')

		data = cls.from_file(input_file, mass_model=mass_model)

		if (input_file.suffix == '.sim' or ''.join(input_file.suffixes) == '.csv.gz') and not output_file.exists():

			data.write_file(output_file)

		elif input_file.suffix == '.hdf5' and not output_file.exists():

			shutil.copy(input_file, output_file)

		return data

	def write_file(self, file):
		'''
		Write ACS data file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to ACS data .hdf5 file
		'''

		acs_data = {'b1': sorted(self.b1, key=lambda x: x[0].to(u.s).value), 'b2': sorted(self.b2, key=lambda x: x[0].to(u.s).value), 'x1': sorted(self.x1, key=lambda x: x[0].to(u.s).value), 
					'x2': sorted(self.x2, key=lambda x: x[0].to(u.s).value), 'y1': sorted(self.y1, key=lambda x: x[0].to(u.s).value), 'y2': sorted(self.y2, key=lambda x: x[0].to(u.s).value)}

		data = {panel: [tuple(value.to(unit).value for value, unit in zip(hit, (u.s, u.keV))) for hit in hits] for panel, hits in acs_data.items()}

		write_hdf5(file, data, file_attributes={'columns': ['time (s)', 'energy (keV)']})

	def plot(self, time_range, energy_range=(80.*u.keV, 2000*u.keV), bin_size=0.05*u.s, file=None, show=False, colors={'b1': 'red', 'b2': 'green', 'x1': 'blue', 'x1': 'orange', 'y1': 'purple', 'y1': 'pink'}, event_time_range=None, title=None):
		'''
		Plot ACS data file.

		Parameters
		----------
		time_range : tuple of astropy.units.Quantity
			Time range of data to plot
		energy_range : tuple of astropy.units.Quantity, optional
			Energy range of data to plot
		bin_size : astropy.units.Quantity, optional
			Time bin size
		file : pathlib.PosixPath, optional
			Path to file to save plot
		show : bool, optional
			Whether to show plot
		colors : dict, optional
			Color to plot for each panel where keys are ACS panel names and entries are colors
		event_time_range : tuple of astropy.units.Quantity, optional
			Time range to highlight on plot
		title : str, optional
			Title of plot
		'''

		nbins = int(np.ceil((time_range[1] - time_range[0]).to_value(u.s) / bin_size.to_value(u.s)))
		duration = nbins * bin_size

		time_range = (time_range[0], time_range[0] + duration)
		bin_edges = np.linspace(time_range[0].to_value(u.s), time_range[1].to_value(u.s), nbins+1)

		for panel in ['b1', 'b2', 'x1', 'x2', 'y1', 'y2']:

			hits = getattr(self, panel)
			times, energies = zip(*hits)

			times = np.array([t.to_value(u.s) for t in times])
			energies = np.array([e.to_value(u.keV) for e in energies])

			mask = (time_range[0].to_value(u.s) <= times) & (times <= time_range[1].to_value(u.s)) & (energy_range[0].to_value(u.keV) <= energies) & (energies <= energy_range[1].to_value(u.keV))

			times = times[mask]
			energies = energies[mask]

			counts, bins = np.histogram(times, bins=bin_edges)
			plt.stairs(counts, bins, color=colors[panel], alpha=0.4)

		plt.legend()

		plt.xlabel('Time (s)')
		plt.ylabel(f'Counts per {bin_size} s bin')

		if event_time_range:
			plt.axvspan(event_time_range[0].to_value(u.s), event_time_range[1].to_value(u.s), facecolor='gray', alpha=0.5)

		if title:
			plt.title(title)

		if file:
			plt.savefig(file)

		if show:
			plt.show()

		plt.close()
