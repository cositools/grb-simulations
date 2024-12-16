import os
import csv
import numpy as np
from .event import event
from .config import read_yaml, define_paths

class source_files():

	def __init__(self, input_file):
		"""
		Create MEGAlib source files from source inputs.

		Parameters
		----------
		input_file : str
			Path to input .yaml file
		"""

		inputs = read_yaml(input_file)

		[self.input_path, self.output_path] = define_paths([inputs['paths']['input'], inputs['paths']['output']], [False, True])

		self.mass_model = inputs['paths']['mass_model']

		if 'shield_counts' in inputs['general'] and inputs['general']['shield_counts'] == 'y' or inputs['general']['shield_counts'] == True:
			self.shield_counts = True
		else:
			self.shield_counts = False

		if 'zenith' in inputs['position']:
			self.zenith = inputs['position']['zenith']
		else:
			self.zenith = None
		if 'zenith_range' in inputs['position']:
			self.zenith_range = inputs['position']['zenith_range']
		else:
			self.zenith_range = None

		if 'azimuth' in inputs['position']:
			self.azimuth = inputs['position']['azimuth']
		else:
			self.azimuth = None
		if 'azimuth_range' in inputs['position']:
			self.azimuth_range = inputs['position']['azimuth_range']
		else:
			self.azimuth_range = None

		if 'latitude' in inputs['position']:
			self.latitude = inputs['position']['latitude']
		else:
			self.latitude = None
		if 'longitude' in inputs['position']:
			self.longitude = inputs['position']['longitude']
		else:
			self.longitude = None

		if 'photon_flux' in inputs['spectra_and_lightcurves']:
			self.ph_flux = inputs['spectra_and_lightcurves']['photon_flux']
		else:
			self.ph_flux = None
		if 'photon_flux_range' in inputs['spectra_and_lightcurves']:
			self.ph_flux_range = inputs['spectra_and_lightcurves']['photon_flux_range']
		else:
			self.ph_flux_range = None

		if 'energy_flux' in inputs['spectra_and_lightcurves']:
			self.e_flux = inputs['spectra_and_lightcurves']['energy_flux']
		else:
			self.e_flux = None
		if 'energy_flux_range' in inputs['spectra_and_lightcurves']:
			self.e_flux_range = inputs['spectra_and_lightcurves']['energy_flux_range']
		else:
			self.e_flux_range = None

		if self.ph_flux == None and self.ph_flux_range == None and self.e_flux == None and self.e_flux_range == None:
			raise RuntimeError("A flux or range of fluxes must be defined in input .yaml file.")

		if 'source_input' in inputs['paths']:
			self.input_source = inputs['paths']['source_input']
		else:
			self.input_source = self.input_path

		if 'orientation' in inputs['paths']:
			self.orientation = inputs['paths']['orientation']
		else:
			self.orientation = None

		if 'source_orientation' in inputs['paths']:
			self.source_orientation = inputs['paths']['source_orientation']
		else:
			self.source_orientation = None
		
		if 'coordinate_system' in inputs['general'] and (inputs['general']['coordinate_system'] == 'local' or inputs['general']['coordinate_system'] == 'galactic'):
			self.coordsys = inputs['general']['coordinate_system']
		else:
			raise RuntimeError("Coordinate system in input .yaml file must be defined as 'local' or 'galactic'.")

		if 'mix_or_match' in inputs['spectra_and_lightcurves']:
			self.mix_or_match = inputs['spectra_and_lightcurves']['mix_or_match']
		else:
			self.mix_or_match = 'match'
		if 'spectrum_type' in inputs['spectra_and_lightcurves']:
			self.spectrum_type = inputs['spectra_and_lightcurves']['spectrum_type']
		else:
			self.spectrum_type = None
		if 'start_time' in inputs['spectra_and_lightcurves']:
			self.start_time = inputs['spectra_and_lightcurves']['start_time']
		elif self.coordsys == 'local':
			self.start_time = None
		else:
			if not self.orientation == None:
				self.orientation_contents = np.loadtxt(self.orientation, usecols=(1, 2, 3, 4, 5), delimiter=' ', skiprows=1, comments=("#", "EN"))
				times = self.orientation_contents[:, 0]
				self.start_time = [times[0] + 50., times[-1] - 50.]
			else:
				raise RuntimeError("An orientation file must be provided to create .source files in galactic coordinates.")

		if self.coordsys == 'galactic':
			if self.latitude == None or self.longitude == None:
				if (self.zenith == None and self.zenith_range == None) or (self.azimuth == None and self.azimuth_range == None):
					raise RuntimeError("Position in local or galactic coordinates must be defined in input .yaml file.")
				else:
					print('Creating .source files in galactic coordinates using source positions provided in local coordinates.')
			elif self.orientation == None:
				raise RuntimeError("Orientation file must be defined in input .yaml file to create .source files in galactic coordinates.")
			else:
				print('Creating .source files in galactic coordinates using source positions provided in galactic coordinates. If local coordinates were provided, these will be ignored.')
		else:
			if (self.zenith == None and self.zenith_range == None) or (self.azimuth == None and self.azimuth_range == None):
				raise RuntimeError("To create .source files in local coordinates, the position must be defined in input .yaml file in local coordinates.")
			else:
				print('Creating .source files in local coordinates.')

	def write_readme(self):
		"""
		Write README for .source file directory.
		"""

		with open(self.output_path + 'README.md', 'w') as f:
			f.write('These source files were created using the input files in ' + self.input_path + '. ')

			if self.coordsys == 'local' or (self.coordsys == 'galactic' and (self.latitude == None or self.longitude == None)):
				if not self.zenith == None:
					if type(self.zenith) == int or type(self.zenith) == float:
						f.write('All .source files have a zenith angle of ' + str(self.zenith) + ' degrees. ')
					elif type(self.zenith) == list:
						f.write('The zenith angles of each source were chosen from the following list (in degrees): ' + str(self.zenith) + '. ')
					else:
						raise RuntimeError("'zenith' in input .yaml file must be int, float, or list.")
				elif not self.zenith_range == None:
					f.write('The zenith angles of each source were chosen randomly between ' + str(self.zenith_range[0]) + ' and ' + str(self.zenith_range[1]) + ' degrees. ')
				else:
					raise RuntimeError("Must specify zenith angle(s) in input .yaml file.")
				if not self.azimuth == None:
					if type(self.azimuth) == int or type(self.azimuth) == float:
						f.write('All .source files have an azimuthal angle of ' + str(azimuth) + ' degrees. ')
					elif type(self.azimuth) == list:
						f.write('The azimuthal angles of each source were chosen from the following list (in degrees): ' + str(self.azimuth) + '. ')
					else:
						raise RuntimeError("'azimuth' in input .yaml file must be int, float, or list.")
				elif not self.azimuth_range == None:
					f.write('The azimuthal angles of each source were chosen randomly between ' + str(self.azimuth_range[0]) + ' and ' + str(self.azimuth_range[1]) + ' degrees. ')
				else:
					raise RuntimeError("Must specify azimuthal angle(s) in input .yaml file.")
			
			else:
				if type(self.latitude) == int or type(self.latitude) == float:
					f.write('All .source files have a galactic latitude of ' + str(self.latitude) + ' degrees. ')
				elif type(self.latitude) == list:
					f.write('The galactic latitudes of each source were chosen randomly between ' + str(self.latitude[0]) + ' and ' + str(self.latitude[1]) + ' degrees. ')
				else:
					raise RuntimeError("'latitude' in input .yaml file must be int, float, or list.")
				if type(self.longitude) == int or type(self.longitude) == float:
					f.write('All .source files have a galactic longitude of ' + str(self.longitude) + ' degrees. ')
				elif type(self.longitude) == list:
					f.write('The galactic longitude of each source were chosen randomly between ' + str(self.longitude[0]) + ' and ' + str(self.longitude[1]) + ' degrees. ')
				else:
					raise RuntimeError("'latitude' in input .yaml file must be int, float, or list.")

			if not self.start_time == None:
				if type(self.start_time) == int or type(self.start_time) == float:
					f.write('All events begin at ' + str(self.start_time) + ' s. ')
				elif type(self.start_time) == list:
					f.write('The start times of each source were chosen randomly between ' + str(self.start_time[0]) + ' and ' + str(self.start_time[1]) + ' s. ')
				else:
					raise RuntimeError("'start_time' in input .yaml file must be int, float, or list.")

			if not self.ph_flux == None:
				if type(self.ph_flux) == int or type(self.ph_flux) == float:
					f.write('All .source files have a flux of ' + str(self.ph_flux) + ' ph/cm<sup>2</sup>/s. ')
				elif type(self.ph_flux) == list:
					f.write('The fluxes of each source were chosen from the following list (in ph/cm<sup>2</sup>/s): ' + str(self.ph_flux) + '. ')
				else:
					raise RuntimeError("'photon_flux' in input .yaml file must be int, float, or list.")
			elif not self.ph_flux_range == None:
				f.write('The fluxes of each source were chosen randomly between ' + str(self.ph_flux_range[0]) + ' and ' + str(self.ph_flux_range[1]) + ' ph/cm<sup>2</sup>/s. ')
			elif not self.e_flux == None:
				if type(self.e_flux) == int or type(self.e_flux) == float:
					f.write('All .source files have a flux of ' + str(e_flux) + ' erg/cm<sup>2</sup>/s. ')
				elif type(self.e_flux) == list:
					f.write('The fluxes of each source were chosen from the following list (in erg/cm<sup>2</sup>/s): ' + str(self.e_flux) + '. ')
				else:
					raise RuntimeError("'energy_flux' in input .yaml file must be int, float, or list.")
			elif not self.e_flux_range == None:
				f.write('The fluxes of each source were chosen randomly between ' + str(self.e_flux_range[0]) + ' and ' + str(self.e_flux_range[1]) + ' erg/cm<sup>2</sup>/s. ')
			else:
				raise RuntimeError("Must specify flux(es) in input .yaml file.")

	def make_event_dict(self):
		"""
		Create list of event names and directories with spectrum & lightcurve file names.
		"""

		key_list = []
		lightcurve_dict = {}
		spectrum_dict = {}
		for file in os.listdir(self.input_path):
			filename = os.fsdecode(file)
			if filename.endswith('_lightcurve.dat'):
				key_list.append(filename[:-15])
				lightcurve_dict[filename[:-15]] = filename

		if self.mix_or_match == 'mix':
			spectrum_list = []
			for file in os.listdir(self.input_path):
				filename = os.fsdecode(file)
				if self.spectrum_type == 'yaml':
					if filename.endswith('_spectrum.yaml'):
						spectrum_list.append(filename)
				elif self.spectrum_type == 'dat':
					if filename.endswith('_spectrum.dat'):
						spectrum_list.append(filename)
				else:
					spectrum_list.append(filename)
			for item in key_list:
				index = np.random.randint(len(spectrum_list))
				filename = spectrum_list[index]
				spectrum_dict[item] = filename
	
		else:
			for item in key_list:
				if self.spectrum_type == 'yaml':
					if os.path.isfile(self.input_path + item + '_spectrum.yaml'):
						spectrum_dict[item] = item + '_spectrum.yaml'
				elif self.spectrum_type == 'dat':
					if os.path.isfile(self.input_path + item + '_spectrum.dat'):
						spectrum_dict[item] = item + '_spectrum.dat'
				else:
					if os.path.isfile(self.input_path + item + '_spectrum.yaml'):
						spectrum_dict[item] = item + '_spectrum.yaml'
					elif os.path.isfile(self.input_path + item + '_spectrum.dat'):
						spectrum_dict[item] = item + '_spectrum.dat'
					elif os.path.isfile(self.input_path + item + '_spectrum.yaml') and os.path.isfile(self.input_path + item + '_spectrum.dat'):
						raise RuntimeError(item + " has both .yaml and .dat spectral files. Spectral type needs to be specified.")
					else:
						raise RuntimeError(item + " has no spectrum file in input directory.")

		self.events = key_list
		self.lightcurves = lightcurve_dict
		self.spectra = spectrum_dict

	def define_spectrum(self, name, spectrum, parameters=None, filename=None, energy_range=None):
		"""
		Determine text to define spectrum in .source file.

		Parameters
		----------
		name : str
			Event name
		spectrum : str
			Spectral type
		parameters : list of int or list of float, optional
			Values of model parameters if spectral type is not 'File'. Length dependent on model
		filename : str, optional
			Path to spectrum if spectral type is 'File'
		energy_range : list of int or list of float, optional
			Low and high energy limits in keV

		Returns
		----------
		spectrum_text : str
			Text to define spectrum in .source file
		"""

		if spectrum == 'File':
			spectrum_text = name + '.Spectrum        File ' + str(filename)
		elif spectrum == 'Mono':
			spectrum_text = name + '.Spectrum        Mono ' + str(parameters[0])
		elif spectrum == 'Band':
			spectrum_text = name + '.Spectrum        Band ' + str(energy_range[0]) + ' ' + str(energy_range[1]) + ' ' + str(parameters[0]) + ' ' + str(parameters[1]) + ' ' + str(parameters[2])
		elif spectrum == 'Comptonized':
			spectrum_text = name + '.Spectrum        Comptonized ' + str(energy_range[0]) + ' ' + str(energy_range[1]) + ' ' + str(parameters[0]) + ' ' + str(parameters[1])
		elif spectrum == 'PowerLaw':
			spectrum_text = name + '.Spectrum        PowerLaw ' + str(energy_range[0]) + ' ' + str(energy_range[1]) + ' ' + str(parameters[0])
		elif spectrum == 'BrokenPowerLaw':
			spectrum_text = name + '.Spectrum        BrokenPowerLaw ' + str(energy_range[0]) + ' ' + str(energy_range[1]) + ' ' + str(parameters[0]) + ' ' + str(parameters[1]) + ' ' + str(parameters[2])
		else:
			raise RuntimeError("Spectral type not supported. 'type' in spectrum .yaml file must be 'Mono', 'Band', 'Comptonized', 'PowerLaw', or 'BrokenPowerLaw'.")

		return spectrum_text

	def lightcurve_spectrum_text(self, name, lightcurve, event_spectrum, source):
		"""
		Read spectra and determine text to define lightcurve & spectrum in .source file.

		Parameters
		----------
		name : str
			Event name
		lightcurve : str
			Lightcurve filename
		event_spectrum : str
			Spectrum filename
		source : grb_simulator.models_sources.event
			Event object for source
		"""

		lightcurve_text = name + '.Lightcurve      File false ' + self.input_source + lightcurve

		if event_spectrum.endswith('_spectrum.yaml'):
			parameters, spectrum, e_range = source.get_spectral_parameters(event_spectrum)
			spectrum_text = self.define_spectrum(name, spectrum['type'], parameters=parameters, energy_range=[spectrum['energy_min'], spectrum['energy_max']])
			return spectrum_text, lightcurve_text, spectrum, parameters, e_range
		elif event_spectrum.endswith('_spectrum.dat'):
			spectrum = {}
			spectrum['type'] = 'File'
			spectrum['filename'] = path + event_spectrum
			spectrum_text = self.define_spectrum(name, spectrum['type'], filename=spectrum['filename'])
			return spectrum_text, lightcurve_text, spectrum
		else:
			raise RuntimeError("Spectrum file name must end in either '_spectrum.yaml' or '_spectrum.dat'.")

	def lightcurve_times(self, lightcurve):
		"""
		Calculate simulation time for .source file based on lightcurve.

		Parameters
		----------
		lightcurve : str
			Lightcurve filename

		Returns
		----------
		sim_time : float
			Maximum time of simulation in s
		"""

		lightcurve_contents = np.loadtxt(self.input_path + lightcurve, usecols=(1, 2), delimiter=' ', skiprows=1, comments=("#", "EN"))

		start_time = np.min(lightcurve_contents[:, 0])
		sim_time = np.max(lightcurve_contents[:, 0]) + 50. - start_time

		return start_time, sim_time

	def edit_lightcurve_times(self, lightcurve, time):
		"""
		Change times in lightcurve file to begin at specified time.

		Parameters
		----------
		lightcurve : str
			Lightcurve filename
		time : float
			Time at which event begins
		"""

		lines = []
		with open(self.input_path + lightcurve, newline='\n') as file:
			reader = csv.reader(file, delimiter=' ', skipinitialspace=True)
			for row in reader:
				lines.append(row)

		lightcurve_contents = np.loadtxt(self.input_path + lightcurve, usecols=(1, 2), delimiter=' ', skiprows=1, comments=("#", "EN"))
		lightcurve_start = np.min(lightcurve_contents[:, 0])
		time_add = time - lightcurve_start

		for i in range(len(lines)):
			if lines[i][0] == 'DP':
				lines[i][1] = str(float(lines[i][1]) + time_add)

		with open(self.input_path + lightcurve, 'w') as file:
			for i in range(len(lines)):
				for j in range(len(lines[i])):
					if j == len(lines[i]) - 1:
						file.write(lines[i][j] + '\n')
					else:
						file.write(lines[i][j] + ' ')

	def make_source_file(self, event, filename, flux, spectrum_text, lightcurve_text, e_flux, sim_time, z=None, a=None, l=None, b=None):
		"""
		Write .source file.

		Parameters
		----------
		event : str
			Event name
		filename : str
			Name of .source file
		flux : int or float
			Flux in ph/cm^2/s
		spectrum_text : str
			Text defining spectrum
		lightcurve_text : str
			Text defining lightcurve
		e_flux : int or float
			Flux in erg/cm^2/s if photon flux was calculated from energy flux
		sim_time : float
			Maximum time of simulation in s
		z : int or float
			Zenith angle in degrees
		a : int or float
			Azimuth angle in degrees
		l : int or float
			Galactic longitude in degrees
		b : int or float
			Galactic latitude in degrees
		"""

		with open(filename, 'w') as f:
			f.write('# Global Parameters\n')
			f.write('Version                     1\n')
			f.write('Geometry                    ' + self.mass_model + '\n')
			f.write('\n# Physics list\n')  
			f.write('PhysicsListEM               LivermorePol\n')
			f.write('\n# Output formats\n')
			f.write('StoreSimulationInfo         init-only\n')
			if self.shield_counts == True:
				f.write('\n# Store shield counts\n')
				f.write('PreTriggerMode              EveryEventWithHits\n')
			f.write('\n# Run and source parameters\n')
			f.write('Run                         GRBSim\n')
			f.write('GRBSim.FileName             ' + event + '\n')
			f.write('GRBSim.Time                 ' + str(sim_time) + '\n')
			if self.coordsys == 'galactic':
				if self.source_orientation == None:
					f.write('GRBSim.OrientationSky       Galactic File NoLoop ' + self.orientation + '\n')
				else:
					f.write('GRBSim.OrientationSky       Galactic File NoLoop ' + self.source_orientation + '\n')
			f.write('GRBSim.Source               ' + event + '\n')
			f.write(event + '.ParticleType    1\n')
			if self.coordsys == 'galactic':
				f.write(event + '.Beam            FarFieldPointSource 0 0 \n')
				if z == None or a == None:
					f.write('\n# Orientation\n')
				else:
					f.write('\n# Orientation corresponding to zenith = ' + str(z) + ' degrees and azimuth = ' + str(a) + ' degrees\n')
				f.write(event + '.Orientation     Galactic Fixed ' + str(b) + ' ' + str(l) + '\n')
			else:
				f.write(event + '.Beam            FarFieldPointSource ' + str(z) + ' ' + str(a) + '\n')
			f.write('\n# Spectrum \n')
			f.write(spectrum_text + '\n')
			if not e_flux == None:
				f.write('\n# Average photon flux in photon/cm2/s corresponding to ' + '{:.2e}'.format(e_flux) + ' erg/cm2/s\n')
			else:
				f.write('\n# Average photon flux in photon/cm2/s\n')
			f.write(event + '.Flux            ' + str(flux) + '\n')
			f.write('\n# Lightcurve\n')
			f.write(lightcurve_text)

	def make_source_files(self):
		"""
		Write .source files for all events in input directory.
		"""

		self.make_event_dict()
		self.write_readme()

		if not self.orientation == None:
			if not hasattr(self, 'orientation_contents'):
				self.orientation_contents = np.loadtxt(self.orientation, usecols=(1, 2, 3, 4, 5), delimiter=' ', skiprows=1, comments=("#", "EN"))
		else:
			self.orientation_contents = None

		for this_event in self.events:
			source_filename = self.output_path + this_event + '.source'
			source = event(self.input_path, self.coordsys, self.orientation_contents)
			source_begin, sim_time = self.lightcurve_times(self.lightcurves[this_event])
			if not self.start_time == None:
				if type(self.start_time) == int or type(self.start_time) == float:
					source_begin = self.start_time
				else:
					source_begin = np.random.uniform(float(self.start_time[0]), float(self.start_time[1]))
				self.edit_lightcurve_times(self.lightcurves[this_event], source_begin)
			if self.spectra[this_event].endswith('_spectrum.yaml'):
				spectrum_text, lightcurve_text, spectrum, parameters, e_range = self.lightcurve_spectrum_text(this_event, self.lightcurves[this_event], self.spectra[this_event], source)
				if self.coordsys == 'local':
					zenith, azimuth, flux, e_flux = source.define_angles_flux(spectrum['type'], parameters, e_range, self.zenith, self.zenith_range, self.azimuth, self.azimuth_range, self.ph_flux, self.ph_flux_range, self.e_flux, self.e_flux_range, self.latitude, self.longitude, source_begin)
					self.make_source_file(this_event, source_filename, flux, spectrum_text, lightcurve_text, e_flux, sim_time+source_begin, z=zenith, a=azimuth)
				else:
					longitude, latitude, flux, e_flux = source.define_angles_flux(spectrum['type'], parameters, e_range, self.zenith, self.zenith_range, self.azimuth, self.azimuth_range, self.ph_flux, self.ph_flux_range, self.e_flux, self.e_flux_range, self.latitude, self.longitude, source_begin)
					self.make_source_file(this_event, source_filename, flux, spectrum_text, lightcurve_text, e_flux, sim_time+source_begin, l=longitude, b=latitude)
			else:
				spectrum_text, lightcurve_text, spectrum = self.lightcurve_spectrum_text(this_event, self.lightcurves[this_event], self.spectra[this_event], source, source_filename)
				raise RuntimeError(".dat spectral files not yet supported.")
