import shutil
import pandas as pd
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import logging
from cosiburstpy.utility.utility import read_hdf5, write_hdf5, SuppressOutput

logger = logging.getLogger(__name__)

class ACSData():

	def __init__(self, data, sort=True, binned=False, bins=None):
		'''
		ACS data handling.

		Parameters
		----------
		data : dict of list of 2-tuple of astropy.units.Quantity
			ACS data where keys are ACS panel names and values are lists of tuples of time and energy if unbinned, and values are 2D np.ndarray if binned
		sort : bool, optional
			Whether to sort by time (ignored if binned)
		binned : bool, optional
			Whether data are binned
		bins : dict of list of astropy.units.Quantity, optional
			Time and energy bin edges
		'''

		self.binned = binned

		if self.binned:

			if bins:

				self.time_bin_edges = bins['time']
				self.energy_bin_edges = bins['energy']

			else:

				raise RuntimeError("Bin edges must be provided for binned data.")

		for panel in data:

			if panel in ['b1', 'b2', 'x1', 'x2', 'y1', 'y2']:

				if sort and not self.binned:
					panel_data = sorted(data[panel], key=lambda x: x[0].to(u.s).value)
				else:
					panel_data = data[panel]

				setattr(self, panel, panel_data)

			else:

				raise RuntimeError(f"{panel} is not a valid panel name.")

	@classmethod
	def from_file(cls, file, mass_model=None, sort=True):
		'''
		Read in ACS data file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to ACS data file
		mass_model : pathlib.PosixPath, optional
			Path to analysis mass model if reading .sim or .sim.gz file
		sort : bool, optional
			Whether to sort by time

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
			raise RuntimeError(f"{file.suffix} files are not supported for ACS data.")

		return acs_data

	@classmethod
	def from_hdf5_file(cls, file, sort=True):
		'''
		Extract ACS hits from .hdf5 file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to ACS data .hdf5 file
		sort : bool, optional
			Whether to sort by time
	
		Returns
		-------
		acs_data : cosiburstpy.acs_data.ACSData
			ACS data
		'''

		data, attributes, dataset_attributes = read_hdf5(file)

		if 'type' in attributes and attributes['type'] == 'binned':

			panels = ['b1', 'b2', 'x1', 'x2', 'y1', 'y2']

			bins = {'time': [value * u.s for value in data['time bins (s)']], 
					'energy': [value * u.keV for value in data['energy bins (keV)']]}
			data = {panel: data[panel] for panel in panels}

			acs_data = cls(data, binned=True, bins=bins)

		else:

			data = {key: [tuple(value * unit for value, unit in zip(hit, (u.s, u.keV))) for hit in hits] for key, hits in data.items()}
			acs_data = cls(data, sort)

		return acs_data

	@classmethod
	def from_sim_file(cls, file, mass_model):
		'''
		Extract ACS hits from .sim or .sim.gz file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to ACS data .sim or .sim.gz file
		mass_model : pathlib.PosixPath
			Path to Data Challenge 3 analysis mass model
	
		Returns
		-------
		acs_data : cosiburstpy.acs_data.ACSData
			ACS data
		'''

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

		logger.info(f"Reading file: {file}")

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
	def combine(cls, files):
		'''
		Combine ACS .hdf5 files.

		Parameters
		----------
		files : list of pathlib.PosixPath
			Paths to binned ACS .hdf5 files

		Returns
		-------
		acs_data : cosiburstpy.acs_data.ACSData
			Combined ACS data
		'''

		logger.info("Combining ACS data files")

		for file in files:
			if file.suffix != '.hdf5':
				raise RuntimeError("Only .hdf5 files can be combined.")

		for i, file in enumerate(files):

			if i == 0:

				data = cls.from_file(file)

				if not data.binned:
					raise RuntimeError("Only binned files can be combined.")

			else:

				component_data = cls.from_file(file)

				if not component_data.binned:
					raise RuntimeError("Only binned files can be combined.")

				if data.time_bin_edges != component_data.time_bin_edges:
					raise RuntimeError("Files must be binned with same time bins.")

				if data.energy_bin_edges != component_data.energy_bin_edges:
					raise RuntimeError("Files must be binned with same energy bins.")

				for panel, value in vars(data).items():
					if panel in ['x1', 'x2', 'y1', 'y2', 'b1', 'b2']:
						setattr(data, panel, np.add(value, getattr(component_data, panel)))

			counts = 0
			for panel, value in vars(data).items():
				if panel in ['x1', 'x2', 'y1', 'y2', 'b1', 'b2']:
					counts += np.sum(value)

			logger.info(f"Added {file}. Total counts: {int(counts)}")

		return data

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

		logger.info(f"Converting {input_file} to .hdf5")

		data = cls.from_file(input_file, mass_model=mass_model)

		if (input_file.suffix == '.sim' or ''.join(input_file.suffixes) == '.csv.gz' or ''.join(input_file.suffixes) == '.sim.gz') and not output_file.exists():

			data.write_file(output_file)

		elif input_file.suffix == '.hdf5' and not output_file.exists():

			shutil.copy(input_file, output_file)

		return data

	def bin(self, time_bins, energy_bins, shift=None):
		'''
		Bin ACS data.

		Parameters
		----------
		time_bins : list of astropy.units.quantity.Quantity
			Time bin edges
		energy_bins : list of astropy.units.quantity.Quantity
			Energy bin edges
		shift : astropy.units.quantity.Quantity, optional
			Time shift to add before binning

		Returns
		-------
		binned_acs_data : cosiburstpy.acs_data.ACSData
			Binned ACS data
		'''

		time_bins = [t.to(u.s).value for t in time_bins]
		energy_bins = [e.to(u.keV).value for e in energy_bins]

		binned_data = {}

		for panel in ['b1', 'b2', 'x1', 'x2', 'y1', 'y2']:

			unbinned_panel_data = getattr(self, panel)
			unbinned_panel_data = [(t.to(u.s).value, e.to(u.keV).value) for t, e in unbinned_panel_data]

			if len(unbinned_panel_data) == 0:

				times = [0.]
				energies = [0.]

			else:

				times, energies = np.array(unbinned_panel_data).T

			if not shift == None:

				for t in range(len(times)):
					times[t] += shift.to(u.s).value

			binned_panel_data, time_edges, energy_edges = np.histogram2d(times, energies, bins=[time_bins, energy_bins])

			binned_data[panel] = binned_panel_data

		binned_acs_data = self.__class__(binned_data, binned=True, bins={'time': time_bins * u.s, 'energy': energy_bins * u.keV})

		return binned_acs_data

	def shift_bins(self, file, time):
		'''
		Update time bins to begin at specified time, and write to file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .hdf5 file
		time : astropy.units.quantity.Quantity
			New start time of data
		'''

		time_add = time - self.time_bin_edges[0]
		self.time_bin_edges = [t + time_add for t in self.time_bin_edges]

		self.write_file(file)

	def write_file(self, file):
		'''
		Write ACS data file or files.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to ACS data .hdf5 file
		'''

		acs_data = {}

		if self.binned:

			dataset_attributes = {}

			for panel in list(vars(self)):

				if panel in ['b1', 'b2', 'x1', 'x2', 'y1', 'y2']:

					acs_data[panel] = getattr(self, panel)

			acs_data['time bins (s)'] = [float(value.to(u.s).value) for value in self.time_bin_edges]
			acs_data['energy bins (keV)'] = [float(value.to(u.keV).value) for value in self.energy_bin_edges]

			write_hdf5(file, acs_data, file_attributes={'type': 'binned'})

		else:

			for panel in list(vars(self)):

				if panel in ['b1', 'b2', 'x1', 'x2', 'y1', 'y2']:

					panel_data = {panel: [tuple(value.to(unit).value for value, unit in zip(hit, (u.s, u.keV))) for hit in getattr(self, panel)]}
					acs_data = acs_data | panel_data

			write_hdf5(file, acs_data, file_attributes={'type': 'unbinned', 'columns': ['time (s)', 'energy (keV)']})

	def plot(self, time_range, energy_range=(80.*u.keV, 2000*u.keV), bin_size=0.05*u.s, file=None, show=False, colors={'b1': 'red', 'b2': 'green', 'x1': 'blue', 'x2': 'orange', 'y1': 'purple', 'y2': 'pink'}, event_time_range=None, title=None, dpi=350):
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
		dpi : int, optional
			Figure resolution
		'''

		nbins = int(np.ceil((time_range[1] - time_range[0]).to_value(u.s) / bin_size.to_value(u.s)))
		duration = nbins * bin_size

		time_range = (time_range[0], time_range[0] + duration)
		bin_edges = np.linspace(time_range[0].to_value(u.s), time_range[1].to_value(u.s), nbins+1)

		for panel in ['b1', 'b2', 'x1', 'x2', 'y1', 'y2']:

			hits = getattr(self, panel)

			if self.binned:

				time_start = np.searchsorted([t.to_value(u.s) for t in self.time_bin_edges], time_range[0].to_value(u.s), side='right') - 1
				time_end = np.searchsorted([t.to_value(u.s) for t in self.time_bin_edges], time_range[1].to_value(u.s), side='right')

				energy_start = np.searchsorted([e.to_value(u.keV) for e in self.energy_bin_edges], energy_range[0].to_value(u.keV), side='right') - 1
				energy_end = np.searchsorted([e.to_value(u.keV) for e in self.energy_bin_edges], energy_range[1].to_value(u.keV) - 0.01, side='right')

				hits = np.array(hits)[time_start:time_end, energy_start:energy_end]
				counts = np.sum(hits, axis=1)

				plt.stairs(counts, [t.to_value(u.s) for t in self.time_bin_edges][time_start:time_end+1], color=colors[panel], alpha=0.4, label=panel)

			else:

				times, energies = zip(*hits)

				times = np.array([t for t in times])
				energies = np.array([e for e in energies])

				mask = (time_range[0].to_value(u.s) <= times) & (times <= time_range[1].to_value(u.s)) & (energy_range[0].to_value(u.keV) <= energies) & (energies <= energy_range[1].to_value(u.keV))

				times = times[mask]
				energies = energies[mask]

				counts, bins = np.histogram(times, bins=bin_edges)
				plt.stairs(counts, bins, color=colors[panel], alpha=0.4, label=panel)

		if self.binned:
			bin_size = round((self.time_bin_edges[1] - self.time_bin_edges[0]).to_value(u.s), 2)

		plt.legend()

		plt.xlabel("Time (s)")
		plt.ylabel(f"Counts per {bin_size} s bin")

		if event_time_range:
			plt.axvspan(event_time_range[0].to_value(u.s), event_time_range[1].to_value(u.s), facecolor='gray', alpha=0.5)

		if title:
			plt.title(title)

		if file:
			plt.savefig(file, dpi=dpi)

		if show:
			plt.show()

		plt.close()
