import numpy as np
import glob
import astropy.units as u
from astropy.stats import bayesian_blocks
import matplotlib.pyplot as plt
from gbm.data import TTE
from gbm.binning.unbinned import bin_by_time
from gbm.background import BackgroundFitter, BackgroundRates
from gbm.background.binned import Polynomial
import gbm.plot
import logging
from cosiburstpy.event.lightcurve import Lightcurve
from cosiburstpy.utility.utility import read_yaml

logger = logging.getLogger(__name__)

class BayesianBlocks():

	def __init__(self, name, data_path, source_time_range, background_time_range, nai_energy_range, bgo_energy_range):
		'''
		Lightcurve binning using Bayesian blocks.

		Parameters
		----------
		name : str
			Name of source
		data_path : pathlib.PosixPath
			Path to directory containing GBM file and TTE data
		source_time_range : 2-tuple of astropy.units.quantity.Quantity
			Time range to consider for source
		background_time_range : 2-tuple of astropy.units.quantity.Quantity
			Time range to fit background
		nai_energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range to include for NaI detectors
		bgo_energy_range : 2-tuple of astropy.units.quantity.Quantity
			Energy range to include for BGO detectors
		'''

		source_data = read_yaml(data_path / f'{name}.yaml')

		self.name = name 

		self.duration = source_data['t90'] * u.s
		self.start_time = source_data['t90_start'] * u.s

		self.nai_paths = glob.glob(str(data_path / '*tte_n*.fit'))
		self.bgo_paths = glob.glob(str(data_path / '*tte_b*.fit'))

		if self.duration <= 0.1 * u.s:
			self.bin_sizes = (0.005 * u.s, 0.002 * u.s)
		else:
			self.bin_sizes = (0.01 * u.s, 0.005 * u.s)

		self.source_time_range = source_time_range
		self.background_time_range = background_time_range

		self.nai_energy_range = nai_energy_range
		self.bgo_energy_range = bgo_energy_range

	def read_tte(self):
		'''
		Read in tte files for each detector and apply time & energy slices for source & background.

		Returns
		-------
		source_tte : gbm.data.phaii.TTE
			TTE data during source interval
		background_tte_nai : gbm.data.phaii.TTE
			TTE data for NaI detectors during background interval
		background_tte_bgo : gbm.data.phaii.TTE
			TTE data for BGO detectors during background interval
		'''

		source_tte = []
		background_tte_nai = []
		background_tte_bgo = []

		for nai_file in self.nai_paths:

			tte = TTE.open(nai_file)

			source_tte_detector = tte.slice_time([self.source_time_range[0].to(u.s).value, self.source_time_range[1].to(u.s).value])
			source_tte_detector = source_tte_detector.slice_energy([self.nai_energy_range[0].to(u.keV).value, self.nai_energy_range[1].to(u.keV).value])
			source_tte.append(source_tte_detector)

			background_tte_detector = tte.slice_time([self.background_time_range[0].to(u.s).value, self.background_time_range[1].to(u.s).value])
			background_tte_nai.append(background_tte_detector)

		for bgo_file in self.bgo_paths:

			tte = TTE.open(bgo_file)

			source_tte_detector = tte.slice_time([self.source_time_range[0].to(u.s).value, self.source_time_range[1].to(u.s).value])
			source_tte_detector = source_tte_detector.slice_energy([self.bgo_energy_range[0].to(u.keV).value, self.bgo_energy_range[1].to(u.keV).value])
			source_tte.append(source_tte_detector)

			background_tte_detector = tte.slice_time([self.background_time_range[0].to(u.s).value, self.background_time_range[1].to(u.s).value])
			background_tte_bgo.append(background_tte_detector)

		source_tte = TTE.merge(source_tte, force_unique=True)
		background_tte_nai = TTE.merge(background_tte_nai, force_unique=True)
		background_tte_bgo = TTE.merge(background_tte_bgo, force_unique=True)

		return source_tte, background_tte_nai, background_tte_bgo

	def bin_tte(self, tte, p0=0.005, spike_tolerance=0.1):
		'''
		Merge and bin data.

		Parameters
		----------
		tte : list of gbm.data.phaii.TTE
			TTE data from each detector
		p0 : float, optional
			False alarm probability to compute prior
		spike_tolerance : float, optional
			Minimum size of first and last bins to include to avoid spikes (s)

		Returns
		-------
		count_rates : np.ndarray
			Count rate in each Bayesian block bin (counts/s)
		bins_bayesian : np.ndarray
			Bayesian block bin edges (s)
		lightcurve_plot : gbm.data.primitives.TimeBins
			Lightcurve for plotting
		nbins : int
			Number of bins for MEGAlib lightcurve
		'''

		times = tte.data.time

		nbins = int((times[-1] - times[0]) / self.bin_sizes[0].to(u.s).value)
		counts, bin_edges = np.histogram(times, bins=nbins)

		bin_centers = []
		for i in range(len(counts)):
			bin_centers.append((bin_edges[i] + bin_edges[i+1]) / 2)

		bins_bayesian = bayesian_blocks(bin_centers, x=counts, p0=p0)
		counts_bayesian, bins_bayesian = np.histogram(times, bins=bins_bayesian)

		if len(counts_bayesian) > 2:

			if bins_bayesian[1] - bins_bayesian[0] < spike_tolerance:

				counts_bayesian = np.delete(counts_bayesian, 0)
				bins_bayesian = np.delete(bins_bayesian, 0)

			if bins_bayesian[-1]-bins_bayesian[-2] < spike_tolerance:

				counts_bayesian = np.delete(counts_bayesian, -1)
				bins_bayesian = np.delete(bins_bayesian, -1)

		bin_sizes = []

		for i in range(len(counts_bayesian)):
			bin_sizes.append(bins_bayesian[i+1] - bins_bayesian[i])

		if min(bin_sizes) > 1.25 * self.duration.to(u.s).value or len(bin_sizes) == 1:
			raise RuntimeError(f"{self.name} not found using Bayesian blocks.")

		count_rates = []

		for i in range(len(counts_bayesian)):
			count_rates.append(counts_bayesian[i] / bin_sizes[i])

		count_rates = np.array(count_rates)

		phaii = tte.to_phaii(bin_by_time, self.bin_sizes[1].to(u.s).value, time_ref=0.0)
		lightcurve_plot = phaii.to_lightcurve()

		return count_rates, bins_bayesian, lightcurve_plot, nbins

	def fit_background(self, background_tte_nai, background_tte_bgo, order=2):
		'''
		Fit and bin background.

		Parameters
		----------
		background_tte_nai : gbm.data.phaii.TTE
			TTE data for NaI detectors during background interval
		background_tte_bgo : gbm.data.phaii.TTE
			TTE data for BGO detectors during background interval
		order : int, optional
			Order of polynomial to fit to background

		Returns
		-------
		background_rates : gbm.background.background.BackgroundRates
			Background rates in NaI and BGO detectors
		'''

		background = []

		binned_nai = background_tte_nai.to_phaii(bin_by_time, self.bin_sizes[0].to(u.s).value, time_ref=0.0)
		binned_bgo = background_tte_bgo.to_phaii(bin_by_time, self.bin_sizes[0].to(u.s).value, time_ref=0.0)

		time_ranges = [(binned_nai.data.tstart[0] + 0.1, self.start_time.to(u.s).value - 1.), (self.duration.to(u.s).value + 1., binned_nai.data.tstop[-1] - 0.1)]

		background_fitter_nai = BackgroundFitter.from_phaii(binned_nai, Polynomial, time_ranges=time_ranges)
		background_fitter_bgo = BackgroundFitter.from_phaii(binned_bgo, Polynomial, time_ranges=time_ranges)

		background_fitter_nai.fit(order=3)
		background_fitter_bgo.fit(order=order)

		background_nai = background_fitter_nai.interpolate_bins(binned_nai.data.tstart, binned_nai.data.tstop)
		background.append(background_nai.integrate_energy(*[self.nai_energy_range[0].to(u.keV).value, self.nai_energy_range[1].to(u.keV).value]))

		background_bgo = background_fitter_bgo.interpolate_bins(binned_nai.data.tstart, binned_nai.data.tstop)
		background.append(background_bgo.integrate_energy(*[self.bgo_energy_range[0].to(u.keV).value, self.bgo_energy_range[1].to(u.keV).value]))

		background_rates = BackgroundRates.sum_time(background)

		return background_rates

	def subtract_background(self, count_rates, bins_bayesian, background_rates):
		'''
		Subtract background from binned lightcurve.

		Parameters
		----------
		count_rates : np.ndarray, shape(N,)
			Count rate in each Bayesian block bin (counts/s)
		bins_bayesian : np.ndarray, shape(N,)
			Bayesian block bin edges (s)
		background_rates : gbm.background.background.BackgroundRates
			Background rates in NaI and BGO detectors

		Returns
		-------
		source_rates : np.ndarray
			Count rate of source in each Bayesian block bin (counts/s)
		'''

		bayesian_rates, bayesian_bins = np.histogram(bins_bayesian[:-1], bins_bayesian, weights=count_rates)

		background_data = []
		for i in range(len(bayesian_rates)):

			index_start = np.inf
			index_end = np.inf

			for j in range(len(background_rates.tstart)):

				if background_rates.tstart[j] <= bayesian_bins[i] and background_rates.tstop[j] > bayesian_bins[i]:

					if background_rates.tstop[j] - bayesian_bins[i] >= bayesian_bins[i] - background_rates.tstart[j]:
						index_start = j
					else:
						index_start = j + 1

				if background_rates.tstart[j] < bayesian_bins[i+1] and background_rates.tstop[j] >= bayesian_bins[i+1]: 
					
					if background_rates.tstop[j] - bayesian_bins[i+1] <= bayesian_bins[i+1] - background_rates.tstart[j]:
						index_end = j + 1
					else:
						index_end = j

			if index_start == np.inf or index_end == np.inf or index_end - index_start == 0:
				raise RuntimeError(f"Background corresponding to Bayesian block bin {i} not found.")

			background_rate = 0
			for j in range(index_start, index_end):
				background_rate += background_rates.rates[j]
			background_rate /= (index_end - index_start)

			background_data.append(background_rate)

		source_rates = bayesian_rates - np.array(background_data)
		source_rates[source_rates < 0] = 0.

		return source_rates

	def plot_lightcurve(self, lightcurve_plot, background_rates, bins_bayesian, count_rates, source_rates, plot_path=None, show=False):
		'''
		Plot lightcurve.

		Parameters
		----------
		lightcurve_plot : gbm.data.primitives.TimeBins
			Lightcurve for plotting
		background_rates : gbm.background.background.BackgroundRates
			Background rates in NaI and BGO detectors
		bins_bayesian : np.ndarray, shape(N,)
			Bayesian block bin edges (s)
		count_rates : np.ndarray, shape(N,)
			Count rate in each Bayesian block bin (counts/s)
		source_rates : np.ndarray, shape(N,)
			Count rate of source in each Bayesian block bin (counts/s)
		plot_path : pathlib.PosixPath, optional
			Path to directory to save plots
		show : bool, optional
			Whether to show plots
		'''

		plot_path.mkdir(parents=True, exist_ok=True)

		lcplot = gbm.plot.Lightcurve(data=lightcurve_plot, background=background_rates)
		lcplot.lightcurve.color='grey'
		lcplot.lightcurve.alpha=0.5
		lcplot.background.color='purple'
		lcplot.background.alpha=0
		lcplot.errorbars.hide()

		plt.axvspan(self.start_time.to(u.s).value, self.duration.to(u.s).value, color='blue', alpha=0.3)

		plt.hist(bins_bayesian[:-1], bins_bayesian, weights=count_rates, color='green', alpha=0.7)
		plt.hist(bins_bayesian[:-1], bins_bayesian, weights=count_rates, histtype='step', color='green', alpha=0.7)

		plt.title(self.name)
		plt.xlim(self.source_time_range[0].to(u.s).value, self.source_time_range[1].to(u.s).value)

		if show:
			plt.show()

		if plot_path:
			plt.savefig(plot_path / f'{self.name}.png')
		
		plt.close()

		bin_centers = []
		bin_widths = []
		for i in range(len(source_rates)):
			bin_centers.append((bins_bayesian[i] + bins_bayesian[i+1]) / 2)
			bin_widths.append(bins_bayesian[i+1] - bins_bayesian[i])

		plt.bar(bin_centers, source_rates, width=bin_widths, color='green', alpha=1)
		plt.axvspan(self.start_time.to(u.s).value, self.duration.to(u.s).value, color='blue', alpha=0.5)

		plt.title(self.name)
		plt.xlim(self.source_time_range[0].to(u.s).value, self.source_time_range[1].to(u.s).value)

		if show:
			plt.show()

		if plot_path:
			plt.savefig(plot_path / f'{self.name}_background_subtracted.png')

		plt.close()

	def rebin_lightcurve(self, bins_bayesian, source_rates, nbins):
		'''
		Rebin lightcurve.

		Parameters
		----------
		bins_bayesian : np.ndarray, shape(N,)
			Bayesian block bin edges (s)
		source_rates : np.ndarray, shape(N,)
			Count rate of source in each Bayesian block bin (counts/s)
		nbins : int
			Number of bins for MEGAlib lightcurve

		Returns
		-------
		lightcurve : cosiburstpy.megalib.lightcurve.Lightcurve
			Lightcurve of source
		'''

		times = list(np.linspace(bins_bayesian[0], bins_bayesian[-1], nbins))
		
		rates = []
		for i in range(len(times) - 1):

			index = np.inf

			for j in range(len(source_rates)):

				if bins_bayesian[j] <= times[i] and bins_bayesian[j+1] > times[i]:

					if bins_bayesian[j+1] - times[i] >= times[i+1] - bins_bayesian[j+1]:
						index = j
					else:
						index = j + 1

			rates.append(source_rates[index])

		first_bin_with_counts = np.inf
		for i in range(len(rates)):
			if first_bin_with_counts > i and rates[i] > 0:
				first_bin_with_counts = i

		times = times[first_bin_with_counts:]
		rates = rates[first_bin_with_counts:]

		times.pop(len(times)-1)

		end = False
		for i, rate in reversed(list(enumerate(rates))):
			if not end:

				if rate == 0:
					times.pop(i)
					rates.pop(i)

				else:
					end = True
					break

		times *= u.s
		rates *= u.s**-1

		start_time = np.min(times)
		time_add = 0. * u.s - start_time

		for i in range(len(times)):
			times[i] += time_add

		lightcurve = Lightcurve(times, rates)

		return lightcurve

	def bin_lightcurve(self, file, plot=False, plot_path=None):
		'''
		Bin lightcurve using Bayesian blocks and create .dat file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .dat file to save lightcurve
		plot_path : pathlib.PosixPath, optional
			Path to directory to save plots
		'''

		file.parent.mkdir(parents=True, exist_ok=True)
		logger.info(f"Creating lightcurve for {self.name}.")

		source_tte, background_tte_nai, background_tte_bgo = self.read_tte()

		try:

			count_rates, bins_bayesian, lightcurve_plot, nbins = self.bin_tte(source_tte)
			background_rates = self.fit_background(background_tte_nai, background_tte_bgo)
			source_rates = self.subtract_background(count_rates, bins_bayesian, background_rates)

			if plot_path:
				self.plot_lightcurve(lightcurve_plot, background_rates, bins_bayesian, count_rates, source_rates, plot_path)

			lightcurve = self.rebin_lightcurve(bins_bayesian, source_rates, nbins)
			lightcurve.write_file(file)

		except:

			logger.warning(f"Error creating lightcurve for {self.name}.")
