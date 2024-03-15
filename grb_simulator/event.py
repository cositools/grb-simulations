import numpy as np
from .config import read_yaml
from .model import model

class event():

	def __init__(self, input_path, coordsys):
		"""
		Define events.

		Parameters
		----------
		input_path : str
			Path to input files
		coordsys : str
			Coordinate system ('local' or 'galactic')
		"""

		self.input_path = input_path
		self.coordsys = coordsys

	def get_spectral_parameters(self, spectrum_file):
		"""
		Read spectral yaml file.

		Parameters
		----------
		spectrum_file : str
			Spectral yaml file

		Returns
		----------
		parameters : list of float
			Spectral parameters
		spectrum : dict
			Contents of spectral yaml file
		e_range : list of float
			Minimum and maximum energies of spectrum in keV
		"""

		parameters = []
		spectrum = read_yaml(self.input_path + spectrum_file)

		if spectrum['type'] == 'Mono':
			parameters.append(spectrum['energy'])
		elif spectrum['type'] == 'Band':
			parameters.append(spectrum['alpha'])
			parameters.append(spectrum['beta'])
			parameters.append(spectrum['ebreak'])
		elif spectrum['type'] == 'Comptonized':
			parameters.append(spectrum['index'])
			parameters.append(spectrum['epeak'])
		elif spectrum['type'] == 'PowerLaw':
			parameters.append(spectrum['index'])
		elif spectrum['type'] == 'BrokenPowerLaw':
			parameters.append(spectrum['ebreak'])
			parameters.append(spectrum['index_lo'])
			parameters.append(spectrum['index_hi'])
			parameters.append(spectrum['energy_max'])
		else:
			raise RuntimeError("Spectral type not supported. 'type' in spectrum .yaml file must be 'Mono', 'Band', 'Comptonized', 'PowerLaw', or 'BrokenPowerLaw'.")

		e_range = [spectrum['energy_min'], spectrum['energy_max']]

		return parameters, spectrum, e_range

	def define_angles_flux(self, spectrum_name=None, parameters=None, e_range=None, zenith=None, zenith_range=None, azimuth=None, azimuth_range=None, ph_flux=None, ph_flux_range=None, e_flux=None, e_flux_range=None):
		"""
		Define zenith angle, azimuth angle, and photon flux of source.

		Parameters
		----------
		spectrum_name : str, optional
			Name of spectral model
		parameters : list of float, optional
			Spectral parameters
		e_range : list of float, optional
			Minimum and maximum energies of spectrum in keV
		zenith : int or float, optional
			Zenith angle
		zenith_range : list of int or list of float, optional
			Minimum and maximum zenith angles
		azimuth : int or float, optional
			Azimuth angle
		azimuth_range : list of int or list of float, optional
			Minimum and maximum azimuth angles
		ph_flux : int or float, optional
			Flux in photons/cm^2/s
		ph_flux_range : list of int or list of float, optional
			Minimum and maximum flux in photons/cm^2/s
		e_flux : int or float, optional
			Flux in ergs/cm^2/s
		e_flux_range : list of int or list of float, optional
			Minimum and maximum flux in ergs/cm^2/s

		Returns
		----------
		zenith : float
			Zenith angle
		azimuth : float
			Azimuth angle
		flux : float
			Photon flux
		e_flux : float
			Energy flux
		"""

		if self.coordsys == 'local':
			if not zenith == None:
				if type(zenith) == int or type(zenith) == float:
					zenith = float(zenith)
				elif type(zenith) == list:
					zenith = float(np.random.choice(zenith))
				else:
					raise RuntimeError("'zenith' in input .yaml file must be int, float, or list.")
			elif not zenith_range == None:
				z_min = float(zenith_range[0])
				z_max = float(zenith_range[1])
				zenith = np.random.uniform(z_min, z_max)
			else:
				raise RuntimeError("Must specify zenith angle(s) in input .yaml file.")

			if not azimuth == None:
				if type(azimuth) == int or type(azimuth) == float:
					azimuth = float(inputs['azimuth'])
				elif type(azimuth) == list:
					azimuth = float(np.random.choice(azimuth))
				else:
					raise RuntimeError("'azimuth' in input .yaml file must be int, float, or list.")
			elif not azimuth_range == None:
				a_min = float(azimuth_range[0])
				a_max = float(azimuth_range[1])
				azimuth = np.random.uniform(a_min, a_max)
			else:
				raise RuntimeError("Must specify azimuthal angle(s) in input .yaml file.")

		else:
			raise RuntimeError("Only detector coordinates are supported for now. 'coordinate_system' in input .yaml file must be 'local'.")

		if not ph_flux == None:
			if type(ph_flux) == int or type(ph_flux) == float:
				flux = float(ph_flux)
			elif type(ph_flux) == list:
				flux = float(np.random.choice(ph_flux))
			else:
				raise RuntimeError("'ph_flux' in input .yaml file must be int, float, or list.")
		elif not ph_flux_range == None:
			flux_min = float(ph_flux_range[0])
			flux_max = float(ph_flux_range[1])
			flux = np.random.uniform(flux_min, flux_max)
		elif not e_flux == None:
			if type(e_flux) == int or type(e_flux) == float:
				pass
			elif type(e_flux) == list:
				e_flux = float(np.random.choice(e_flux))
			else:
				raise RuntimeError("'e_flux' in input .yaml file must be int, float, or list.")
		elif not e_flux_range == None:
			e_flux_min = float(e_flux_range[0])
			e_flux_max = float(e_flux_range[1])
			e_flux = np.random.uniform(e_flux_min, e_flux_max)
		else:
			raise RuntimeError("Must specify flux(es) in input .yaml file.")

		if not e_flux == None:
			spectral_model = model(spectrum_name)
			flux = float(spectral_model.calc_photon_flux(e_flux, e_range, parameters))

		self.zenith = zenith
		self.azimuth = azimuth
		self.flux = flux
		self.e_flux = e_flux

		return zenith, azimuth, flux, e_flux


