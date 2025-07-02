from scipy import integrate
import astropy.units as u
from cosiburstpy.utility.utility import read_yaml, write_yaml
from .spectral_models import mono, band, comp, pl, bpl, sbpl

class Spectrum():

	def __init__(self, model, parameters):
		'''
		Define spectral model.

		Parameters
		----------
		model : str
			Name of model (either 'Mono', 'Band', 'Comptonized', 'PowerLaw', 'BrokenPowerLaw', or 'SmoothlyBrokenPowerLaw')
		parameters : dict of astropy.units.quantity.Quantity
			Values of model parameters where keys are parameter names and values are parameter values
		'''

		self.model = model 

		if parameters:
			self.set_model(parameters)

	def set_model(self, parameters):
		'''
		Set spectral model.

		Parameters
		----------
		parameters : dict of astropy.units.quantity.Quantity, optional
			Values of model parameters where keys are parameter names and values are parameter values
		'''

		self.parameters = parameters

		if self.model == 'Mono':

			self.parameter_list = [float(f"{self.parameters['energy'].to(u.keV).value:.2f}")]
			self.function = lambda e: mono(e, self.parameters['energy'])

		elif self.model == 'Band':

			if not 'ebreak' in self.parameters:
				self.parameters['ebreak'] = parameters['epeak'] / (parameters['alpha'] + 2)

			if not 'epeak' in self.parameters:
				self.parameters['epeak'] = self.parameters['ebreak'] * (parameters['alpha'] + 2)

			self.parameter_list = [float(f"{self.parameters['alpha']:.2f}"), float(f"{self.parameters['beta']:.2f}"), float(f"{self.parameters['epeak'].to(u.keV).value:.2f}")]
			self.function = lambda e: band(e, self.parameters['alpha'], self.parameters['beta'], self.parameters['ebreak'])

		elif self.model == 'Comptonized':

			self.parameter_list = [float(f"{self.parameters['index']:.2f}"), float(f"{self.parameters['epeak'].to(u.keV).value:.2f}")]
			self.function = lambda e: comp(e, self.parameters['index'], self.parameters['epeak'])

		elif self.model == 'PowerLaw':

			self.parameter_list = [float(f"{self.parameters['index']:.2f}")]
			self.function = lambda e: pl(e, self.parameters['index'])

		elif self.model == 'BrokenPowerLaw':

			self.parameter_list = [float(f"{self.parameters['ebreak'].to(u.keV).value:.2f}"), float(f"{self.parameters['index_lo']:.2f}"), float(f"{self.parameters['index_hi']:.2f}"), float(f"{self.parameters['emax'].to(u.keV).value:.2f}")]
			self.function = lambda e: bpl(e, self.parameters['ebreak'], self.parameters['index_lo'], self.parameters['index_hi'], self.parameters['emax'])

		elif self.model == 'SmoothlyBrokenPowerLaw':

			self.parameter_list = [float(f"{self.parameters['ebreak'].to(u.keV).value:.2f}"), float(f"{self.parameters['index_lo']:.2f}"), float(f"{self.parameters['index_hi']:.2f}"), float(f"{self.parameters['bscale']:.2f}")]
			self.function = lambda e: sbpl(e, self.parameters['ebreak'], self.parameters['index_lo'], self.parameters['index_hi'], self.parameters['bscale'])

		else:

			raise RuntimeError("Spectral model must be 'Mono', 'PowerLaw', 'BrokenPowerLaw', 'Band', 'Comptonized', or 'SmoothlyBrokenPowerLaw'.")

	@classmethod
	def from_file(cls, file):
		'''
		Read in spectrum file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .yaml file

		Returns
		-------
		spectrum : cosiburstpy.event.spectrum.Spectrum
			Spectrum
		'''

		data = read_yaml(file)
		spectral_model = data.pop('type')

		for key in data:
			if key in ['energy', 'epeak', 'ebreak', 'emax']:
				data[key] *= u.keV

		spectrum = Spectrum(spectral_model, data)

		return spectrum

	def integrate_photon_spectrum(self, energy_range):
		'''
		Integrate photon spectrum over energy range.

		Parameters
		----------
		energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range

		Returns
		-------
		integral : astropy.units.quantity.Quantity
			Value of integral
		'''

		if self.model == 'Mono':

			if energy_range[0] <= self.parameters['energy'] <= energy_range[1]:
				integral = 1. * u.dimensionless_unscaled
			else:
				integral = 0. * u.dimensionless_unscaled

		else:

			function = lambda e: self.function(e * u.keV).to(u.keV**-1).value
			integral = integrate.quad(function, energy_range[0].to(u.keV).value, energy_range[1].to(u.keV).value)[0] * u.dimensionless_unscaled

		return integral

	def integrate_energy_spectrum(self, energy_range):
		'''
		Integrate energy spectrum over energy range.

		Parameters
		----------
		energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range

		Returns
		-------
		integral : astropy.units.quantity.Quantity
			Value of integral
		'''

		if self.model == 'Mono':

			if energy_range[0] <= self.parameters['energy'] <= energy_range[1]:
				integral = self.parameters['energy'].to(u.erg)
			else:
				integral = 0. * u.erg

		else:

			function = lambda e: (e * u.keV * self.function(e * u.keV)).to(u.dimensionless_unscaled).value
			integral = integrate.quad(function, energy_range[0].to(u.keV).value, energy_range[1].to(u.keV).value)[0] * u.keV

		return integral

	def photon_flux_from_energy_flux(self, energy_flux, energy_range):
		'''
		Calculate photon flux from energy flux.

		Parameters
		----------
		energy_flux : astropy.units.quantity.Quantity
			Energy flux
		energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range

		Returns
		-------
		photon_flux : astropy.units.quantity.Quantity
			Photon flux
		'''

		photon_integral = self.integrate_photon_spectrum(energy_range)
		energy_integral = self.integrate_energy_spectrum(energy_range)

		if not energy_integral == 0:

			photon_flux = (energy_flux * photon_integral / energy_integral).to(u.cm**-2 * u.s**-1)

		elif energy_flux == 0.:

			photon_flux = 0. * u.cm**-2 * u.s**-1

		else:

			raise RuntimeError("Integral of energy spectrum is 0.")

		return photon_flux

	def energy_flux_from_photon_flux(self, photon_flux, energy_range):
		'''
		Calculate energy flux from photon flux.

		Parameters
		----------
		photon_flux : astropy.units.quantity.Quantity
			Photon flux
		energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range

		Returns
		-------
		energy_flux : astropy.units.quantity.Quantity
			Energy flux
		'''

		photon_integral = self.integrate_photon_spectrum(energy_range)
		energy_integral = self.integrate_energy_spectrum(energy_range)

		if not photon_integral == 0:

			energy_flux = (photon_flux * energy_integral / photon_integral).to(u.erg / u.cm**2 / u.s)

		elif photon_flux == 0.:

			energy_flux = 0. * u.erg / u.cm**2 / u.s

		else:

			raise RuntimeError("Integral of photon spectrum is 0.")

		return energy_flux

	def change_photon_flux_energy_range(self, old_photon_flux, old_energy_range, new_energy_range):
		'''
		Calculate photon flux in new energy range.

		Parameters
		----------
		old_photon_flux : astropy.units.quantity.Quantity
			Photon flux
		old_energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range of input flux
		new_energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range of output flux

		Returns
		-------
		new_photon_flux : astropy.units.quantity.Quantity
			Photon flux in new energy range
		'''

		old_photon_integral = self.integrate_photon_spectrum(old_energy_range)
		new_photon_integral = self.integrate_photon_spectrum(new_energy_range)

		new_photon_flux = old_photon_flux * new_photon_integral / old_photon_integral

		return new_photon_flux

	def change_energy_flux_energy_range(self, old_energy_flux, old_energy_range, new_energy_range):
		'''
		Calculate energy flux in new energy range.

		Parameters
		----------
		old_energy_flux : astropy.units.quantity.Quantity
			Energy flux
		old_energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range of input flux
		new_energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range of output flux

		Returns
		-------
		new_energy_flux : astropy.units.quantity.Quantity
			Energy flux in new energy range
		'''

		old_energy_integral = self.integrate_energy_spectrum(old_energy_range)
		new_energy_integral = self.integrate_energy_spectrum(new_energy_range)

		new_energy_flux = old_energy_flux * new_energy_integral / old_energy_integral

		return new_energy_flux

	def generate_text(self, energy_range=None):
		'''
		Generate text to define spectrum in MEGAlib .source file.

		Parameters
		----------
		energy_range : 2-tuple of astropy.units.quantity.Quantity, optional
			Energy range of simulation

		Returns
		-------
		spectrum_text : str
			Text to define spectrum in .source file
		'''

		if not energy_range and self.model != 'Mono':
			raise RuntimeError("Energy range must be provided if spectral model is not 'Mono'.")

		if self.model == 'Mono':

			spectrum_text = f"Mono {self.parameters['energy'].to(u.keV).value:.2f}"

		elif self.model == 'Band':

			spectrum_text = f"Band {energy_range[0].to(u.keV).value:.2f} {energy_range[1].to(u.keV).value:.2f} {self.parameters['alpha']:.2f} {self.parameters['beta']:.2f} {self.parameters['ebreak'].to(u.keV).value:.2f}"

		elif self.model == 'Comptonized':

			spectrum_text = f"Comptonized {energy_range[0].to(u.keV).value:.2f} {energy_range[1].to(u.keV).value:.2f} {self.parameters['index']:.2f} {self.parameters['epeak'].to(u.keV).value:.2f}"

		elif self.model == 'PowerLaw':

			spectrum_text = f"PowerLaw {energy_range[0].to(u.keV).value:.2f} {energy_range[1].to(u.keV).value:.2f} {self.parameters['index']:.2f}"

		elif self.model == 'BrokenPowerLaw':

			spectrum_text = f"BrokenPowerLaw {energy_range[0].to(u.keV).value:.2f} {energy_range[1].to(u.keV).value:.2f} {self.parameters['ebreak'].to(u.keV).value:.2f} {self.parameters['index_lo']:.2f} {self.parameters['index_hi']:.2f}"

		else:

			raise RuntimeError("Spectral model must be 'Mono', 'PowerLaw', 'BrokenPowerLaw', 'Band', or 'Comptonized'.")

		return spectrum_text

	def write_file(self, file):
		'''
		Write spectrum to .yaml file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .yaml file to save spectrum
		'''

		file.parent.mkdir(parents=True, exist_ok=True)

		spectrum = {}
		spectrum['type'] = self.model

		for key in self.parameters:

			if isinstance(self.parameters[key], u.quantity.Quantity):
				spectrum[key] = float(self.parameters[key].value)
			else:
				spectrum[key] = self.parameters[key]

		write_yaml(file, spectrum)
