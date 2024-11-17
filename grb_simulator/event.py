import numpy as np
from astropy.coordinates import SkyCoord, spherical_to_cartesian
import astropy.units as u
from scoords import Attitude, SpacecraftFrame
from .config import read_yaml
from .model import model

class event():

	def __init__(self, input_path, coordsys, orientation=None):
		"""
		Define events.

		Parameters
		----------
		input_path : str
			Path to input files
		coordsys : str
			Coordinate system ('local' or 'galactic')
		orientation : np.array, optional
			Contents of orientation file
		"""

		self.input_path = input_path
		self.coordsys = coordsys
		self.orientation = orientation

		if not type(self.orientation) == 'NoneType':
			self.t = self.orientation[:, 0]
			self.xb = self.orientation[:, 1]
			self.xl = self.orientation[:, 2]
			self.zb = self.orientation[:, 3]
			self.zl = self.orientation[:, 4]

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

	def calculate_galactic_position(self, zenith, azimuth, x_latitude, x_longitude, z_latitude, z_longitude):
		"""
		Calculate position in galactic coordinates from spacecraft pointing, zenith, and azimuth.

		Parameters
		----------
		zenith : float
			Zenith angle in degrees
		azimuth : float
			Azimuth angle in degrees
		x_latitude : float
			Galactic latitude in degrees of spacecraft x direction
		x_longitude : float
			Galactic longitude in degrees of spacecraft x direction
		z_latitude : float
			Galactic latitude in degrees of spacecraft z direction
		z_longitude : float
			Galactic longitude in degrees of spacecraft z direction

		Returns
		----------
		longitude : float
			Galactic longitude of source in degrees
		latitude : float
			Galactic latitude of source in degrees
		"""

		x = SkyCoord(np.deg2rad(x_longitude)*u.rad, np.deg2rad(x_latitude)*u.rad, frame='galactic')
		z = SkyCoord(np.deg2rad(z_longitude)*u.rad, np.deg2rad(z_latitude)*u.rad, frame='galactic')
		attitude = Attitude.from_axes(x=x, z=z, frame='galactic')
		
		x_loc, y_loc, z_loc = spherical_to_cartesian(1, np.deg2rad(zenith)*u.rad, np.deg2rad(azimuth)*u.rad)
		c_loc = SkyCoord(x_loc, y_loc, z_loc, representation_type='cartesian', frame = SpacecraftFrame(attitude=attitude))
		c_galactic = c_loc.transform_to('galactic')
		
		longitude = c_galactic.l.deg
		latitude = c_galactic.b.deg	

		return longitude, latitude

	def define_angles_flux(self, spectrum_name=None, parameters=None, e_range=None, zenith=None, zenith_range=None, azimuth=None, azimuth_range=None, ph_flux=None, ph_flux_range=None, e_flux=None, e_flux_range=None, latitude=None, longitude=None, time=None):
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
			Zenith angle in degrees
		zenith_range : list of int or list of float, optional
			Minimum and maximum zenith angles
		azimuth : int or float, optional
			Azimuth angle in degrees
		azimuth_range : list of int or list of float, optional
			Minimum and maximum azimuth angles in degrees
		ph_flux : int or float, optional
			Flux in photons/cm^2/s
		ph_flux_range : list of int or list of float, optional
			Minimum and maximum flux in photons/cm^2/s
		e_flux : int or float, optional
			Flux in ergs/cm^2/s
		e_flux_range : list of int or list of float, optional
			Minimum and maximum flux in ergs/cm^2/s
		latitude : int, float, or list, optional
			Galactic latitude or minimum and maximum galactic latitudes in degrees
		longitude : int, float, or list, optional
			Galactic longitude or minimum and maximum galactic longitudes in degrees
		time : float, optional
			Time at which event begins if calculating galactic position from zenith and azimuth

		Returns
		----------
		zenith or longitude : float
			Zenith angle (for local coordinate system) or galactic longitude (for galactic coordinate system) in degrees
		azimuth or latitude : float
			Azimuth angle (for local coordinate system) or galactic latitude (for galactic coordinate system) in degrees
		flux : float
			Photon flux in photons/cm^2/s
		e_flux : float
			Energy flux in ergs/cm^2/s
		"""

		if self.coordsys == 'local' or latitude == None or longitude == None:
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
					azimuth = float(azimuth)
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

			if self.coordsys == 'galactic':
				index = np.absolute(self.t - time).argmin()
				longitude, latitude = self.calculate_galactic_position(zenith, azimuth, self.xb[index], self.xl[index], self.zb[index], self.zl[index])

		else:
			if type(latitude) == float or type(latitude) == int:
				latitude = float(latitude)
			elif type(latitude) == list:
				latitude = np.random.uniform(float(latitude[0]), float(latitude[1]))

			if type(longitude) == float or type(longitude) == int:
				longitude = float(longitude)
			elif type(longitude) == list:
				longitude = np.random.uniform(float(longitude[0]), float(longitude[1]))
			
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

		self.flux = flux
		self.e_flux = e_flux

		if self.coordsys == 'local':
			self.zenith = zenith
			self.azimuth = azimuth
			return zenith, azimuth, flux, e_flux
		else:
			self.longitude = longitude
			self.latitude = latitude
			return longitude, latitude, flux, e_flux
