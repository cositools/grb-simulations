import numpy as np
import astropy.units as u
from .position import check_if_earth_occulted

class Event():

	def __init__(self, name, lightcurve, spectrum, orientation=None):
		'''
		Define events.

		Parameters
		----------
		name : str
			Event name
		lightcurve : cosiburstpy.Lightcurve
			Event lightcurve
		spectrum : cosiburstpy.event.spectrum.Spectrum
			Spectral model
		orientation : cosiburstpy.simulations.spacecraft_orientation.SpacecraftOrientation, optional
			Spacecraft orientation during event
		'''

		self.name = name

		self.lightcurve = lightcurve
		self.time_range = (np.min(self.lightcurve.bin_edges), np.max(self.lightcurve.bin_edges))
		self.duration = self.time_range[1] - self.time_range[0]

		self.spectrum = spectrum

		if orientation:
			self.orientation = orientation

			if hasattr(orientation, 'exclude'):
				self.excluded = any(orientation.exclude)

	def set_average_flux(self, flux, energy_range):
		'''
		Set average flux.

		Parameters
		----------
		flux : astropy.units.quantity.Quantity
			Average photon or energy flux
		energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range of flux
		'''

		self._energy_range = energy_range

		if flux.unit.is_equivalent(u.cm**-2 * u.s**-1):

			self._average_flux = flux

		elif flux.unit.is_equivalent(u.erg / u.cm**2 / u.s):

			self._average_flux = self.spectrum.photon_flux_from_energy_flux(flux, energy_range)

		else:

			raise RuntimeError(f"{flux.unit.to_string()} is an invalid flux unit. Must be equivalent to 1/cm2/s for photon flux or erg/cm2/s for energy flux.")

		self._peak_flux = self.peak_flux_from_average_flux(self._average_flux)

	def set_peak_flux(self, flux, energy_range):
		'''
		Set peak flux.

		Parameters
		----------
		flux : astropy.units.quantity.Quantity
			Peak photon or energy flux
		energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range of flux
		'''

		self._energy_range = energy_range

		if flux.unit.is_equivalent(u.cm**-2 * u.s**-1):

			self._peak_flux = flux

		elif flux.unit.is_equivalent(u.erg / u.cm**2 / u.s):

			self._peak_flux = self.spectrum.photon_flux_from_energy_flux(flux, energy_range)

		else:

			raise RuntimeError(f"{flux.unit.to_string()} is an invalid flux unit. Must be equivalent to 1/cm2/s for photon flux or erg/cm2/s for energy flux.")

		self._average_flux = self.average_flux_from_peak_flux(self._peak_flux)

	@property
	def energy_range(self):
		'''
		Energy range.
		'''

		if not hasattr(self, '_energy_range'):
			raise RuntimeError("No energy range defined. Must set using set_average_flux() or set_peak_flux().")

		return self._energy_range

	def average_photon_flux(self, energy_range=None):
		'''
		Average photon flux.

		Parameters
		----------
		energy_range : 2-tuple of astropy.units.quantity.Quantity, optional
			Energy range of flux

		Returns
		-------
		average_flux : astropy.units.quantity.Quantity
			Average photon flux
		'''

		if not hasattr(self, '_average_flux'):	
			raise RuntimeError("No flux defined. Must set using set_average_flux() or set_peak_flux().")

		if not energy_range:
			average_flux = self._average_flux
		else:
			average_flux = self.spectrum.change_photon_flux_energy_range(self._average_flux, self._energy_range, energy_range)
			
		return average_flux

	def peak_photon_flux(self, energy_range=None):
		'''
		Peak photon flux.

		Parameters
		----------
		energy_range : 2-tuple of astropy.units.quantity.Quantity, optional
			Energy range of flux

		Returns
		-------
		peak_flux : astropy.units.quantity.Quantity
			Peak photon flux
		'''

		if not hasattr(self, '_peak_flux'):
			raise RuntimeError("No flux defined. Must set using set_average_flux() or set_peak_flux()")

		if not energy_range:
			peak_flux = self._peak_flux
		else:
			peak_flux = self.spectrum.change_photon_flux_energy_range(self._peak_flux, self._energy_range, energy_range)

		return peak_flux

	def average_energy_flux(self, energy_range=None):
		'''
		Average energy flux.

		Parameters
		----------
		energy_range : 2-tuple of astropy.units.quantity.Quantity, optional
			Energy range of flux

		Returns
		-------
		average_flux : astropy.units.quantity.Quantity
			Average energy flux
		'''

		if not hasattr(self, '_average_flux'):
			raise RuntimeError("No flux defined. Must set using set_average_flux() or set_peak_flux()")

		energy_flux = self.spectrum.energy_flux_from_photon_flux(self._average_flux, self._energy_range)

		if not energy_range:
			average_flux = energy_flux
		else:
			average_flux = self.spectrum.change_energy_flux_energy_range(energy_flux, self._energy_range, energy_range)

		return average_flux

	def peak_energy_flux(self, energy_range=None):
		'''
		Peak energy flux.

		Parameters
		----------
		energy_range : 2tuple of astropy.units.quantity.Quantity, optional
			Energy range of flux

		Returns
		-------
		peak_flux : astropy.units.quantity.Quantity
			Peak energy flux
		'''

		if not hasattr(self, '_peak_flux'):
			raise RuntimeError("No flux defined. Must set using set_average_flux() or set_peak_flux()")

		energy_flux = self.spectrum.energy_flux_from_photon_flux(self._peak_flux, self._energy_range)

		if not energy_range:
			peak_flux = energy_flux
		else:
			peak_flux = self.spectrum.change_energy_flux_energy_range(energy_flux, self._energy_range, energy_range)

		return peak_flux

	def peak_flux_from_average_flux(self, flux):
		'''
		Calculate peak flux from average flux.

		Parameters
		----------
		flux : astropy.units.quantity.Quantity
			Average photon or energy flux

		Returns
		-------
		peak_flux : astropy.units.quantity.Quantity
			Peak photon or energy flux
		'''

		peak_flux = flux * self.lightcurve.peak_ratio

		return peak_flux

	def average_flux_from_peak_flux(self, flux):
		'''
		Calculate average flux from peak flux.

		Parameters
		----------
		flux : astropy.units.quantity.Quantity
			Peak photon or energy flux

		Returns
		-------
		average_flux : astropy.units.quantity.Quantity
			Average photon or energy flux
		'''

		average_flux = flux / self.lightcurve.peak_ratio

		return average_flux

	@property
	def position(self):
		'''
		Source position.
		'''

		if not hasattr(self, '_position'):
			raise RuntimeError("Position not yet defined.")

		return self._position 

	@position.setter
	def position(self, position):
		'''
		Set position.

		Parameters
		----------
		position : astropy.coordinates.sky_coordinate.SkyCoord
			Source position
		'''

		self._position = position

		if hasattr(self, 'orientation'):

			occultation_start = check_if_earth_occulted(self.position, self.orientation.earth_zeniths[0], self.orientation.altitudes[0])
			occultation_stop = check_if_earth_occulted(self.position, self.orientation.earth_zeniths[-1], self.orientation.altitudes[-1])

			if occultation_start or occultation_stop:
				self.occulted = True
			else:
				self.occulted = False

	@property
	def polarization(self):
		'''
		Polarization fraction and angle.
		'''

		if not hasattr(self, '_position'):
			raise RuntimeError("Polarization not yet defined.")

		return self._polarization

	@polarization.setter
	def polarization(self, polarization):
		'''
		Set polarization.

		Parameters
		----------
		polarization : 2-tuple of float and astropy.coordinates.angles.core.Angle
			Polarization in the form (polarization fraction, polarization angle in IAU convention)
		'''

		self._polarization = polarization

