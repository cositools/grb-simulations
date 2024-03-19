import numpy as np
from scipy import integrate

class model():

	def __init__(self, model):
		"""
		Define spectral models.

		Parameters
		----------
		model : str
			Name of model (for now, either 'Mono', 'Band', 'Comptonized', 'PowerLaw', or 'BrokenPowerLaw')
		"""

		self.model = model 

	def band(self, e, alpha, beta, ebreak):
		"""
		Band function.

		Parameters
		----------
		e : float
			Energy (keV)
		alpha : float
			Low energy spectral index
		beta : float
			High energy spectral index
		ebreak : float
			Break energy

		Returns
		----------
		amplitude : float
			Amplitude of spectrum at given energy
		"""

		if e <= (alpha - beta) * ebreak:
			amplitude = (e / 100)**alpha * np.exp(-e / ebreak)
		else:
			amplitude = (e / 100)**beta * np.exp(beta - alpha) * ((alpha - beta) * ebreak / 100)**(alpha - beta)
		
		return amplitude

	def comp(self, e, index, epeak):
		"""
		Comptonized function.

		Parameters
		----------
		e : float
			Energy (keV)
		index : float
			Spectral index
		epeak : float
			Peak energy

		Returns
		----------
		amplitude : float
			Amplitude of spectrum at given energy
		"""

		amplitude = e**index * np.exp(-(index + 2) * e / epeak)

		return amplitude

	def pl(self, e, index):
		"""
		Power law function.

		Parameters
		----------
		e : float
			Energy (keV)
		index : float
			Spectral index

		Returns
		----------
		amplitude : float
			Amplitude of spectrum at given energy
		"""

		amplitude = e**(-index)

		return amplitude

	def bpl(self, e, ebreak, index_lo, index_hi, e_max):
		"""
		Broken power law function.

		Parameters
		----------
		e : float
			Energy (keV)
		ebreak : float
			Break energy (keV)
		index_lo : float
			Low energy spectral index
		index_hi : float
			High energy spectral index
		emax : float
			Maximum energy (keV)

		Returns
		----------
		amplitude : float
			Amplitude of spectrum at given energy
		"""

		if e <= ebreak:
			amplitude = e**(-index_lo)
		else:
			amplitude = e**(-index_hi) * emax**(index_hi - index_lo)

		return amplitude

	def set_model(self, parameters):
		"""
		Set spectral model.

		Parameters
		----------
		parameters : list
			Values of model parameters. Length dependent on model
		"""

		self.parameters = parameters

		if self.model == 'Mono':
			self.parameter_names = ['energy']
			self.function = lambda e: e
		elif self.model == 'Band':
			self.parameter_names = ['alpha', 'beta', 'break energy']
			self.function = lambda e: self.band(e, parameters[0], parameters[1], parameters[2])
		elif self.model == 'Comptonized':
			self.parameter_names = ['index', 'peak energy']
			self.function = lambda e: self.comp(e, parameters[0], parameters[1])
		elif self.model == 'PowerLaw':
			self.parameter_names = ['index']
			self.function = lambda e: self.pl(e, parameters[0])
		elif self.model == 'BrokenPowerLaw':
			self.parameter_names = ['break energy', 'low index', 'high index', 'max energy']
			self.function = lambda e: self.bpl(e, parameters[0], parameters[1], parameters[2], parameters[3])
		else:
			raise RuntimeError("Spectral model not supported. Must be 'Mono', 'Band', 'Comptonized', 'PowerLaw', or 'BrokenPowerLaw'.")

		if len(self.parameters) != len(self.parameter_names):
			raise RuntimeError(f"The length of the parameters list does not match the number of model parameters. "
							   f"The parameters for the " + self.model + " model are " + str(self.parameter_names) + '.')

	def calc_photon_flux(self, e_flux, e_range=None, parameters=None):
		"""
		Calculate photon flux from energy flux.

		Parameters
		----------
		e_flux : float
			Energy flux in erg/s/cm^2
		e_range : list of int or list of float, optional
			Low and high energy limits in keV
		parameters : list of int or list of float, optional
			Values of model parameters. Length dependent on model

		Returns
		----------
		ph_flux : float
			Photon flux in photons/s/cm^2
		"""

		erg_to_kev = 6.24150647e8
		e_times_function = lambda e: e * self.function(e)

		if not parameters == None:
			self.set_model(parameters)

		if self.model == 'Mono':
			ph_flux = e_flux * erg_to_kev / self.parameters[0]
		else:
			if not e_range == None:
				ph_flux = e_flux * erg_to_kev * integrate.quad(self.function, e_range[0], e_range[1])[0] / integrate.quad(e_times_function, e_range[0], e_range[1])[0]
			else:
				raise RuntimeError("Must provide energy range for integration.")

		return ph_flux

	def calc_energy_flux(self, ph_flux, e_range=None, parameters=None):
		"""
		Calculate energy flux from photon flux.

		Parameters
		----------
		ph_flux : float
			Photon flux in photons/s/cm^2
		e_range : list of int or list of float, optional
			Low and high energy limits in keV
		parameters : list of int or list of float, optional
			Values of model parameters. Length dependent on model

		Returns
		----------
		e_flux : float
			Energy flux in erg/s/cm^2
		"""

		erg_to_kev = 6.24150647e8
		e_times_function = lambda e: e * self.function(e)

		if not parameters == None:
			self.set_model(parameters)

		if self.model == 'Mono':
			e_flux = ph_flux * self.parameters[0] / erg_to_kev
		else:
			if not e_range == None:
				e_flux = ph_flux * integrate.quad(e_times_function, e_range[0], e_range[1])[0] / integrate.quad(self.function, e_range[0], e_range[1])[0] / erg_to_kev
			else:
				raise RuntimeError("Must provide energy range for integration.")

		return e_flux