import astropy.units as u
import logging
from cosiburstpy.utility.utility import read_yaml, write_yaml
from cosiburstpy.event.spectrum import Spectrum

logger = logging.getLogger(__name__)

class DefineSpectrum():

	def __init__(self, name, gbm_file):
		'''
		Define spectral model for source.

		Parameters
		----------
		name : str
			Name of source
		gbm_file : pathlib.PosixPath
			Path to .yaml file containing GBM burst catalog data for source
		'''

		self.name = name

		self.source_data = read_yaml(gbm_file)

		if self.source_data['flnc_best_fitting_model'] == 'flnc_band':

			parameters = {'alpha': self.source_data['flnc_band_alpha'], 'beta': self.source_data['flnc_band_beta'], 'ebreak': self.source_data['flnc_band_epeak'] / (self.source_data['flnc_band_alpha'] + 2) * u.keV}
			self.spectrum = Spectrum('Band', parameters)

		elif self.source_data['flnc_best_fitting_model'] == 'flnc_comp':

			parameters = {'index': self.source_data['flnc_comp_index'], 'epeak': self.source_data['flnc_comp_epeak'] * u.keV}
			self.spectrum = Spectrum('Comptonized', parameters)

		elif self.source_data['flnc_best_fitting_model'] == 'flnc_plaw':

			parameters = {'index': self.source_data['flnc_plaw_index']}
			self.spectrum = Spectrum('PowerLaw', parameters)

		elif self.source_data['flnc_best_fitting_model'] == 'flnc_sbpl':

			parameters = {'index_lo': self.source_data['flnc_sbpl_indx1'], 'index_hi': self.source_data['flnc_sbpl_indx2'], 'ebreak': self.source_data['flnc_sbpl_brken'] * u.keV, 'bscale': self.source_data['flnc_sbpl_brksc']}
			self.spectrum = Spectrum('SmoothlyBrokenPowerLaw', parameters)

		else:

			logger.warning(f'No spectral fit available for {self.name}.')
			self.spectrum = None

	def save_spectrum(self, file):
		'''
		Create spectral .yaml file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .yaml file to save spectrum
		'''

		logger.info(f'Creating spectrum for {self.name}.')

		self.spectrum.write_file(file)
