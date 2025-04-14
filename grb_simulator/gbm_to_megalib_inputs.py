import numpy as np
import glob
import matplotlib.pyplot as plt
import os
import astropy.stats
from gbm.data import TTE
from gbm.binning.unbinned import bin_by_time
from gbm.background import BackgroundFitter, BackgroundRates
from gbm.background.binned import Polynomial
from gbm.plot import Lightcurve
from .config import read_yaml, write_yaml, define_paths

class gbm_to_megalib_inputs():

	def __init__(self, input_file):
		"""
		Convert data downloaded from GBM to MEGAlib input files.

		Parameters
		----------
		input_file : str
			Path to input .yaml file
		"""

		inputs = read_yaml(input_file)

		[self.input_path, 
		 self.output_path, 
		 self.plot_path] = define_paths([inputs['paths']['input'], inputs['paths']['output'], inputs['paths']['plots']], 
		 								[False, True, True])

		self.energy_range = inputs['energy']['cosi_range']
		self.background_t_range = inputs['time']['background_range']
		self.source_t_range = inputs['time']['source_range']
		self.nai_e_range = inputs['energy']['nai_range']
		self.bgo_e_range = inputs['energy']['bgo_range']

		self.n_sources = self.count_sources()

	def count_sources(self):
		"""
		Count number of events in input directory.

		Returns
		----------
		n : int
			Number of events
		"""

		n = 0
		for name in os.listdir(self.input_path):
			if os.path.isdir(self.input_path + name):
				n += 1

		return n

	def lightcurve_binning(self, t90):
		"""
		Determine time bin size based on event duration.

		Parameters
		----------
		t90 : float
			Duration (T90) of event

		Returns
		----------
		bin_rate : float
			Time bin size for plotting
		bin_rate2 : float
			Time bin size
		"""

		if t90 <= 0.1:
			bin_rate = 0.002
			bin_rate2 = 0.005
		elif t90 > 0.1 and t90 <= 1.0:
			bin_rate = 0.005
			bin_rate2 = 0.01
		else:
			bin_rate = 0.005
			bin_rate2 = 0.01

		return bin_rate, bin_rate2

	def read_detectors(self, detect, merge, energy_range):
		"""
		Read in detector files and apply time & energy slices.

		Parameters
		----------
		detect : list of str
			Detector files
		merge : list of gbm.data.phaii.TTE
			TTE data from each detector
		energy_range : list of float
			Minimum and maximum energies for detector

		Returns
		----------
		bkgd_list : list
			TTE data from each detector to use for background fitting
		merge : list
			TTE data from each detector
		"""

		bkgd_list = []               # empty array to hold/manipulate background data

		for i in range(len(detect)):
			tte = TTE.open(detect[i]) # data
			bkgd = tte.slice_time(self.background_t_range) # background
			tte = tte.slice_time(self.source_t_range) # data time slice
			tte = tte.slice_energy(energy_range) # data energy slice
			bkgd_list.append(bkgd)
			merge.append(tte) # merge data

		return bkgd_list, merge

	def bin_data(self, merge, bin_rate, bin_rate2, t90):
		"""
		Merge and bin data.

		Parameters
		----------
		merge : list of gbm.data.phaii.TTE
			TTE data from each detector
		bin_rate : float
			Time bin size for plotting and background fitting
		bin_rate2 : float
			Time bin size for MEGAlib lightcurve
		t90 : float
			Duration (T90) of event

		Returns
		----------
		rates : list of float
			Count rate in each bin
		bins2 : np.ndarray
			Bayesian block bins
		lc_data : gbm.data.primitives.TimeBins
			Lightcurve for plotting
		nbins : int
			Number of bins for MEGAlib lightcurve
		err : bool
			Whether there was an issue finding event using Bayesian blocks
		"""

		tte_merge = TTE.merge(merge, force_unique=True) # merge NaI and BGO data
		t = tte_merge.data.time # assign time values
		nbins = int((t[-1] - t[0]) / bin_rate2) # number of bins for initial binning
		counts, bins, __ = plt.hist(t, bins=nbins) # bin data initially
		bin_centers = []
		for j in range(len(counts)):
			bin_centers.append((bins[j+1]+bins[j])/2)
		bins2 = astropy.stats.bayesian_blocks(bin_centers, x=counts, p0=0.005) # bin using Bayesian blocks
		counts2, bins2, __ = plt.hist(t, bins=bins2)
		plt.clf()
		if len(counts2) > 2: # get rid of spiky first & last bins
			if bins2[1]-bins2[0] < 0.1:
				counts2 = np.delete(counts2, 0)
				bins2 = np.delete(bins2, 0)
			if bins2[-1]-bins2[-2] < 0.1:
				counts2 = np.delete(counts2, -1)
				bins2 = np.delete(bins2, -1)
		rates = []
		for j in range(len(counts2)):
			rates.append(counts2[j]/(bins2[j+1]-bins2[j]))
		phaii = tte_merge.to_phaii(bin_by_time, bin_rate, time_ref=0.0) # create time binned phaii data for plotting
		lc_data = phaii.to_lightcurve() # create data lightcurve for plotting

		bin_sizes = np.array([])
		for j in range(len(counts2)):
			bin_sizes = np.append(bin_sizes, bins2[j+1] - bins2[j])

		err = False
		if min(bin_sizes) > t90 or len(bin_sizes) == 1:
			print('Not finding event with Bayesian blocks. Not saving')
			err = True

		return rates, bins2, lc_data, nbins, err

	def fit_background(self, n_detect, b_detect, bkgd_list_nai, bkgd_list_bgo, bin_rate, t90_start, t90):
		"""
		Fit and bin background.

		Parameters
		----------
		n_detect : list of str
			NaI detector files
		b_detect : list of str
			BGO detector files
		bkgd_list_nai : list of gbm.data.phaii.TTE
			TTE data from each NaI detector
		bkgd_list_bgo : list of gbm.data.phaii.TTE
			TTE data from each BGO detector
		bin_rate : float
			Time bin size
		t90_start : float
			Start time of event
		t90 : float
			Duration of event

		Returns
		----------
		tot_bkgd : list of gbm.background.background.BackgroundRates
			Background rates in NaI and BGO detectors
		err : bool
			Whether there was an issue fitting background
		"""
		
		tot_bkgd = []                               # empty array for summing NaI/BGO background data

		# Merge and bin background
		if len(n_detect) > 0:
			bkgd_merge_nai = TTE.merge(bkgd_list_nai, force_unique=True) # merge NaI background data
			bkgd_phaii_nai = bkgd_merge_nai.to_phaii(bin_by_time, bin_rate, time_ref=0.0) # create time binned phaii data for NaI background
		if len(b_detect) > 0:
			bkgd_merge_bgo = TTE.merge(bkgd_list_bgo, force_unique=True) # merge BGO background data
			bkgd_phaii_bgo = bkgd_merge_bgo.to_phaii(bin_by_time, bin_rate, time_ref=0.0) # create time binned phaii data for BGO background

		# Perform background fit, interpolate, integrate over energy, and sum
		src_range = [(bkgd_phaii_nai.data.tstart[0]+0.1, t90_start-0.5), (t90+1.0, bkgd_phaii_nai.data.tstop[-1]-0.1)] # NaI background time range

		err = False
		try:
			backfitter_n = BackgroundFitter.from_phaii(bkgd_phaii_nai, Polynomial, time_ranges=src_range) # initialize
			backfitter_n.fit(order=2) # NaI background fit
			bkgd_n = backfitter_n.interpolate_bins(bkgd_phaii_nai.data.tstart, bkgd_phaii_nai.data.tstop) # NaI backround interpolation
			bkgd_n = bkgd_n.integrate_energy(*self.nai_e_range) # NaI backround energy integration
			tot_bkgd.append(bkgd_n)
			backfitter_b = BackgroundFitter.from_phaii(bkgd_phaii_bgo, Polynomial, time_ranges=src_range) # initialize
			backfitter_b.fit(order=2) # BGO background fit
			bkgd_b = backfitter_b.interpolate_bins(bkgd_phaii_nai.data.tstart, bkgd_phaii_nai.data.tstop) # BGO backround interpolation
			bkgd_b = bkgd_b.integrate_energy(*self.bgo_e_range) # BGo backround energy integration
			tot_bkgd.append(bkgd_b)
		except:
			print('Error fitting background. Not saving')
			err = True

		return tot_bkgd, err

	def subtract_background(self, bb_data, lc_bkgd):
		"""
		Subtract background from lightcurve.

		Parameters
		----------
		bb_data : tuple
			Lightcurve binned using Bayesian blocks
		lc_bkgd : gbm.background.background.BackgroundRates
			Background rates

		Returns
		----------
		data_centers : list of float
			Centers of time bins
		data_heights : list of float
			Count rates in bins
		data_widths : list of float
			Widths of time bins
		"""

		# Rebin background into same bins as data
		bkgd_data = []
		for j in range(len(bb_data[0])):
			index_start = np.inf
			index_end = np.inf
			for k in range(len(lc_bkgd.tstart)):
				if lc_bkgd.tstart[k] <= bb_data[1][j] and lc_bkgd.tstop[k] > bb_data[1][j]: # find first background bin that overlaps with data bin
					if lc_bkgd.tstop[k] - bb_data[1][j] >= bb_data[1][j] - lc_bkgd.tstart[k]: # most of background bin is within data bin
						index_start = k
					else:
						index_start = k+1
				if lc_bkgd.tstart[k] < bb_data[1][j+1] and lc_bkgd.tstop[k] >= bb_data[1][j+1]: # find last background bin that overlaps with data bin
					if lc_bkgd.tstop[k] - bb_data[1][j+1] <= bb_data[1][j+1] - lc_bkgd.tstart[k]: # most of background bin is within data bin
						index_end = k+1
					else:
						index_end = k
			tot_rate = 0
			if index_end - index_start == 0:
				avg_rate = lc_bkgd.rates[index_start]
			elif index_start != np.inf and index_end != np.inf: 
				for k in range(index_start, index_end):
					tot_rate += lc_bkgd.rates[k]
				avg_rate = tot_rate / (index_end - index_start) 
			bkgd_data.append(avg_rate)

		# Subtract background
		source_rate = bb_data[0] - bkgd_data
		for j in range(len(source_rate)):
			if source_rate[j] < 0.0:
				source_rate[j] = 0.0

		data_centers = []
		data_heights = []
		data_widths = []
		for j in range(len(source_rate)):
			start = float(bb_data[1][j])
			end = float(bb_data[1][j+1])
			data_centers.append((start+end)/2)
			data_heights.append(float(source_rate[j]))
			data_widths.append(end-start)

		return data_centers, data_heights, data_widths

	# Plot lightcurve
	def plot_lightcurve(self, lc_data, lc_bkgd, bins, rates, data_centers, data_heights, data_widths, t90_start, t90, name):
		"""
		Subtract background from lightcurve.

		Parameters
		----------
		lc_data : gbm.data.primitives.TimeBins
			Source+background lightcurve
		lc_bkgd : gbm.background.background.BackgroundRates
			Background rates
		bins : np.ndarray
			Bayesian block bins
		rates : list of float
			Count rate in each bin
		data_centers : list of float
			Centers of time bins
		data_heights : list of float
			Count rates in bins
		data_widths : list of float
			Widths of time bins
		t90_start : float
			Start time of event
		t90 : float
			Duration of event
		name : str
			Name of event
		"""

		# Plot the lightcurves with the selections and background fit and overlay Bayesian Blocks
		lcplot = Lightcurve(data=lc_data, background=lc_bkgd)
		lcplot.lightcurve.color='grey'
		lcplot.lightcurve.alpha=0.5
		lcplot.background.color='purple'
		lcplot.background.alpha=0
		lcplot.errorbars.hide()

		plt.axvspan(t90_start, t90, color='blue', alpha=0.3)

		plt.hist(bins[:-1], bins, weights=rates, color='green', alpha=0.7)
		plt.hist(bins[:-1], bins, weights=rates, histtype='step', color='green', alpha=0.7)

		plt.title(name)
		plt.xlim(self.source_t_range[0], self.source_t_range[1])
		plt.savefig(self.plot_path + '%s.png' % name)
		plt.clf()

		# Plot background subtracted lightcurve
		plt.bar(data_centers, data_heights, width=data_widths, color='green', alpha=1)
		plt.axvspan(t90_start, t90, color='blue', alpha=0.5)
		plt.title(name)
		plt.xlim(self.source_t_range[0], self.source_t_range[1])
		plt.savefig(self.plot_path + '%s_bkg_subtracted.png' % name)
		plt.clf()

	def write_lightcurve(self, name, bb_data, nbins, data_heights):
		"""
		Write lightcurve to file.

		Parameters
		----------
		bb_data : tuple
			Lightcurve binned using Bayesian blocks
		nbins : int
			Number of bins for MEGAlib lightcurve
		data_heights : list of float
			Count rates in bins
		"""

		start_list = []
		rate_list = []

		lc_bins = np.linspace(bb_data[1][0], bb_data[1][-1], nbins)

		begin = False
		for j in range(len(lc_bins) - 1):
			index = np.inf
			for k in range(len(bb_data[0])):
				if bb_data[1][k] <= lc_bins[j] and bb_data[1][k+1] > lc_bins[j]: # find first BB bin that overlaps with smaller bin
					if bb_data[1][k+1] - lc_bins[j] >= lc_bins[j+1] - lc_bins[j]:
						index = k
					else:
						index = k+1
						index_start = k
			rate = data_heights[index]

			if not begin:
				if not float(rate) == 0:
					begin = True
					start_list.append(float(lc_bins[j]))
					rate_list.append(float(rate))
			else:
				start_list.append(float(lc_bins[j]))
				rate_list.append(float(rate))

		end = False
		for i, e in reversed(list(enumerate(rate_list))):
			if not end:
				if e == 0:
					start_list.pop(i)
					rate_list.pop(i)
				else:
					end = True
					break

		with open(self.output_path + name + '_lightcurve.dat', 'w') as f:
			f.write('# DP | Time (s) | Count rate (counts/s) \n')
			f.write('# IP LinLin\n')
			for i in range(len(start_list)):
				f.write('DP ' + str(start_list[i] - start_list[0]) + ' ' + str(rate_list[i]) + '\n')
			f.write('EN')

	def make_lightcurve(self, name, source_info):
		"""
		Create lightcurve.

		Parameters
		----------
		name : str
			Name of event
		source_info : dict
			GBM data for event

		Returns
		----------
		err : bool
			Whether there was an issue fitting background
		"""

		print('Making lightcurve...')

		t_centers = []

		bin_rate, bin_rate2 = self.lightcurve_binning(float(source_info['t90']))

		# Find all detector files in event folder
		n_detect = glob.glob(self.input_path + name + '/*tte_n*.fit')
		b_detect = glob.glob(self.input_path + name + '/*tte_b*.fit')

		# Reset background data
		bkgd_phaii_nai = 0
		bkgd_phaii_bgo = 0

		merge = []                          # empty array to merge detector data
		bkgd_list_nai, merge = self.read_detectors(n_detect, merge, self.nai_e_range)
		bkgd_list_bgo, merge = self.read_detectors(b_detect, merge, self.bgo_e_range)

		rates, bins2, lc_data, nbins, err = self.bin_data(merge, bin_rate, bin_rate2, float(source_info['t90']))
		if not err:
			tot_bkgd, err = self.fit_background(n_detect, b_detect, bkgd_list_nai, bkgd_list_bgo, bin_rate, float(source_info['t90_start']), float(source_info['t90']))

		if not err:
			lc_bkgd = BackgroundRates.sum_time(tot_bkgd) # sum NaI and BGO backgrounds
			bb_data = plt.hist(bins2[:-1], bins2, weights=rates, color='green', alpha=0.7)

			data_centers, data_heights, data_widths = self.subtract_background(bb_data, lc_bkgd)

			self.plot_lightcurve(lc_data, lc_bkgd, bins2, rates, data_centers, data_heights, data_widths, source_info['t90_start'], source_info['t90'], source_info['trigger_name'])
			self.write_lightcurve(source_info['trigger_name'], bb_data, nbins, data_heights)

		return err

	def make_spectrum(self, source_info):
		"""
		Create spectrum.

		Parameters
		----------
		source_info : dict
			GBM data for event

		Returns
		----------
		err : bool
			Whether best fit spectrum is unsupported
		"""

		print('Making spectrum file...')

		err = False

		if source_info['flnc_best_fitting_model'] == 'flnc_band':
			band_dict = {}
			band_dict['type'] = 'Band'
			band_dict['energy_min'] = self.energy_range[0]
			band_dict['energy_max'] = self.energy_range[1]
			band_dict['alpha'] = source_info['flnc_band_alpha']
			band_dict['beta'] = source_info['flnc_band_beta']
			band_dict['ebreak'] = source_info['flnc_band_epeak'] / (source_info['flnc_band_alpha'] + 2) 
			write_yaml(self.output_path + source_info['trigger_name'] + '_spectrum.yaml', band_dict)
		elif source_info['flnc_best_fitting_model'] == 'flnc_comp':
			comptonized_dict = {}
			comptonized_dict['type'] = 'Comptonized'
			comptonized_dict['energy_min'] = self.energy_range[0]
			comptonized_dict['energy_max'] = self.energy_range[1]
			comptonized_dict['index'] = source_info['flnc_comp_index']
			comptonized_dict['epeak'] = source_info['flnc_comp_epeak']
			write_yaml(self.output_path + source_info['trigger_name'] + '_spectrum.yaml', comptonized_dict)
		else:
			print('For now, only using events with Band function or Comptonized best spectral fits. Not saving')
			err = True

		return err

	def make_spectrum_lightcurve(self, name):
		"""
		Make spectrum and lightcurve for event.

		Parameters
		----------
		name : str
			Name of event
		"""

		source_info = read_yaml(self.input_path + name + '/' + name + '.yaml')
		if (os.path.isfile(self.output_path + source_info['trigger_name'] + '_spectrum.yaml') and 
			os.path.isfile(self.output_path + source_info['trigger_name'] + '_lightcurve.dat')):
				print('Spectrum and lightcurve already exist')
		else:
			if os.path.isfile(self.output_path + source_info['trigger_name'] + '_spectrum.yaml'):
				print('Spectrum file exists already. Deleting old file')
				os.remove(self.output_path + source_info['trigger_name'] + '_spectrum.yaml')
			elif os.path.isfile(self.output_path + source_info['trigger_name'] + '_lightcurve.dat'):
				print('Lightcurve file exists already. Deleting old file')
				os.remove(self.output_path + source_info['trigger_name'] + '_lightcurve.dat')
			s_error = self.make_spectrum(source_info)
			if not s_error:
				lc_error = self.make_lightcurve(name, source_info)
				if lc_error:
					print('Deleting spectrum')
					os.remove(self.output_path + source_info['trigger_name'] + '_spectrum.yaml')

	def write_readme(self):
		"""
		Write README for MEGAlib input file directory.
		"""

		with open(self.output_path + 'README.md', 'w') as f:
			f.write('# Spectra and Lightcurves\n\n')
			f.write('This directory contains spectra and lightcurves of the Fermi-GBM data found in `' + self.input_path + 
					'`. The lightcurves and spectra were binned using Bayesian blocks. The energy range for the spectral .yaml files is ' + 
					str(self.energy_range[0]) + '-' + str(self.energy_range[1]) + ' keV. Only GRBs with Band and Comptonized best spectral fits are included.')

	def make_spectra_lightcurves(self):
		"""
		Make spectra and lightcurves for all events in directory.
		"""

		self.write_readme()
		
		n = 0
		for name in os.listdir(self.input_path):
			if os.path.isdir(self.input_path + name):
				n += 1
				print(str(n) + '/' + str(self.n_sources) + ': ' + name)
				self.make_spectrum_lightcurve(name)
