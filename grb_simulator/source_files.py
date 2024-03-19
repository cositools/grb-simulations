import os
import numpy as np
from .event import event
from .config import define_paths

class source_files():

	def __init__(self, inputs):
		"""
		Create MEGAlib source files from source inputs.

		Parameters
		----------
		inputs : dict
			Contents of input yaml file 
		input_path : str
			Path to MEGAlib input files
		output_path : str
			Path to save MEGAlib source files
		"""

		[self.input_path, self.output_path] = define_paths([inputs['input_path'], inputs['output_path']], [False, True])

		if 'shield_counts' in inputs.keys() and inputs['shield_counts'] == 'y' or inputs['shield_counts'] == True:
			self.shield_counts = True
		else:
			self.shield_counts = False

		if 'zenith' in inputs.keys():
			self.zenith = inputs['zenith']
		else:
			self.zenith = None
		if 'zenith_min' in inputs.keys() and 'zenith_max' in inputs.keys():
			self.zenith_range = [inputs['zenith_min'], inputs['zenith_max']]
		else:
			self.zenith_range = None
		if self.zenith == None and self.zenith_range == None:
			raise RuntimeError("A zenith angle or range of angles must be defined in input .yaml file.")

		if 'azimuth' in inputs.keys():
			self.azimuth = inputs['azimuth']
		else:
			self.azimuth = None
		if 'azimuth_min' in inputs.keys() and 'azimuth_max' in inputs.keys():
			self.azimuth_range = [inputs['azimuth_min'], inputs['azimuth_max']]
		else:
			self.azimuth_range = None
		if self.azimuth == None and self.azimuth_range == None:
			raise RuntimeError("An azimuth angle or range of angles must be defined in input .yaml file.")

		if 'ph_flux' in inputs.keys():
			self.ph_flux = inputs['ph_flux']
		else:
			self.ph_flux = None
		if 'ph_flux_min' in inputs.keys() and 'ph_flux_max' in inputs.keys():
			self.ph_flux_range = [inputs['ph_flux_min'], inputs['ph_flux_max']]
		else:
			self.ph_flux_range = None

		if 'e_flux' in inputs.keys():
			self.e_flux = inputs['e_flux']
		else:
			self.e_flux = None
		if 'e_flux_min' in inputs.keys() and 'e_flux_max' in inputs.keys():
			self.e_flux_range = [inputs['e_flux_min'], inputs['e_flux_max']]
		else:
			self.e_flux_range = None

		if self.ph_flux == None and self.ph_flux_range == None and self.e_flux == None and self.e_flux_range == None:
			raise RuntimeError("A flux or range of fluxes must be defined in input .yaml file.")

		if 'source_input_path' in inputs.keys():
			self.input_source = inputs['source_input_path']
		else:
			self.input_source = None
		
		if 'mass_model_path' in inputs.keys():
			self.mass_model = inputs['mass_model_path']
		else:
			raise RuntimeError("Mass model path must be defined in input .yaml file.")
		if inputs['coordinate_system'] == 'local':
			self.coordsys = inputs['coordinate_system']
		else:
			raise RuntimeError("Only detector coordinates are supported for now. 'coordinate_system' in input .yaml file must be 'local'.")

		if 'mix_or_match' in inputs.keys():
			self.mix_or_match = inputs['mix_or_match']
		else:
			self.mix_or_match = 'match'
		if 'spectrum_type' in inputs.keys():
			self.spectrum_type = inputs['spectrum_type']
		else:
			self.spectrum_type = None

	def write_readme(self):
		"""
		Write README for source file directory.
		"""

		with open(self.output_path + 'README.md', 'w') as f:
			f.write('These source files were created using the input files in ' + self.input_path + '. ')
		
			if self.coordsys == 'local':
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
				raise RuntimeError("Only detector coordinates are supported for now. 'coordinate_system' in input .yaml file must be 'local'.")

			if not self.ph_flux == None:
				if type(self.ph_flux) == int or type(self.ph_flux) == float:
					f.write('All .source files have a flux of ' + str(self.ph_flux) + ' ph/cm<sup>2</sup>/s. ')
				elif type(self.ph_flux) == list:
					f.write('The fluxes of each source were chosen from the following list (in ph/cm<sup>2</sup>/s): ' + str(self.ph_flux) + '. ')
				else:
					raise RuntimeError("'ph_flux' in input .yaml file must be int, float, or list.")
			elif not self.ph_flux_range == None:
				f.write('The fluxes of each source were chosen randomly between ' + str(self.ph_flux_range[0]) + ' and ' + str(self.ph_flux_range[1]) + ' ph/cm<sup>2</sup>/s. ')
			elif not self.e_flux == None:
				if type(self.e_flux) == int or type(self.e_flux) == float:
					f.write('All .source files have a flux of ' + str(e_flux) + ' erg/cm<sup>2</sup>/s. ')
				elif type(self.e_flux) == list:
					f.write('The fluxes of each source were chosen from the following list (in erg/cm<sup>2</sup>/s): ' + str(self.e_flux) + '. ')
				else:
					raise RuntimeError("'e_flux' in input .yaml file must be int, float, or list.")
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
		Determine text to define spectrum in source file.

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
			Text to define spectrum in source file
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
		Read spectra and determine text to define lightcurve & spectrum in source file.

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

		if self.input_source == None:
			path = self.input_path
		else:
			path = self.input_source

		lightcurve_text = name + '.Lightcurve      File false ' + path + lightcurve

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

	def make_source_file(self, event, filename, z, a, flux, spectrum_text, lightcurve_text, e_flux):
		"""
		Write source file.

		Parameters
		----------
		event : str
			Event name
		filename : str
			Soure file
		z : int or float
			Zenith angle
		a : int or float
			Azimuth angle
		flux : int or float
			Flux in ph/cm^2/s
		spectrum_text : str
			Text defining spectrum
		lightcurve_text : str
			Text defining lightcurve
		e_flux : int or float
			Flux in erg/cm^2/s if photon flux was calculated from energy flux
		"""

		with open(filename, 'w') as f:
			f.write('# Global Parameters\n')
			f.write('Version                     1\n')
			f.write('Geometry                    ' + self.mass_model + '\n')
			f.write('\n# Physics list\n')  
			f.write('PhysicsListEM               LivermorePol\n')
			f.write('\n# Output formats\n')
			f.write('StoreSimulationInfo         init-only\n')
			if self.shield_counts == False:
				f.write('\n# Store shield counts\n')
				f.write('PreTriggerMode              EveryEventWithHits\n')
			f.write('\n# Run and source parameters\n')
			f.write('Run                         GRBSim\n')
			f.write('GRBSim.FileName             ' + event + '\n')
			f.write('GRBSim.Time                 5000.0\n')
			f.write('GRBSim.Source               ' + event + '\n')
			f.write(event + '.ParticleType    1\n')
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
		Write source files for all events in input directory.
		"""

		self.make_event_dict()

		for this_event in self.events:
			source_filename = self.output_path + this_event + '.source'
			source = event(self.input_path, self.coordsys)
			if self.spectra[this_event].endswith('_spectrum.yaml'):
				spectrum_text, lightcurve_text, spectrum, parameters, e_range = self.lightcurve_spectrum_text(this_event, self.lightcurves[this_event], self.spectra[this_event], source)
				zenith, azimuth, flux, e_flux = source.define_angles_flux(spectrum['type'], parameters, e_range, self.zenith, self.zenith_range, self.azimuth, self.azimuth_range, self.ph_flux, self.ph_flux_range, self.e_flux, self.e_flux_range)
				self.make_source_file(this_event, source_filename, zenith, azimuth, flux, spectrum_text, lightcurve_text, e_flux)
			else:
				spectrum_text, lightcurve_text, spectrum = self.lightcurve_spectrum_text(this_event, self.lightcurves[this_event], self.spectra[this_event], source, source_filename)
				raise RuntimeError(".dat spectral files not yet supported.")
