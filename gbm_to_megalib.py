import numpy as np
import glob
from gbm.data import TTE
from gbm.binning.unbinned import bin_by_time
from gbm.background import BackgroundFitter
from gbm.background.binned import Polynomial
from gbm.background import BackgroundRates
from gbm.plot import Lightcurve
from gbm.plot import Spectrum
import gbm.plot
import astropy.stats
import matplotlib.pyplot as plt
import csv
import yaml
import os
import argparse

# Initialize parser & read arguments from command line
def parse_args():

	parser = argparse.ArgumentParser()
	parser.add_argument("-y", "--yaml", help = "Path to input .yaml file", default='examples/example_gbm_to_megalib.yaml')
	args = parser.parse_args()

	return args

# Read yaml file
def read_yaml(file):

	with open(file, "r") as myfile:
		inputs = yaml.safe_load(myfile)

	return inputs

# Define path to output files
def define_paths(inputs):

	input_path = os.path.abspath(inputs['input_path']) + '/'
	output_path = os.path.abspath(inputs['output_path']) + '/'
	if not os.path.isdir(output_path):
		os.makedirs(output_path)
	plot_path = os.path.abspath(inputs['plot_path']) + '/'
	if not os.path.isdir(plot_path):
		os.makedirs(plot_path)

	return input_path, output_path, plot_path

# Write yaml file with data for each burst (update to save spectra)
def write_yaml(inputs, data, output_path):

	mydict = {}
	for i in range(len(inputs['download'])):
		if type(data[i]) == np.float64:
			mydict[inputs['download'][i]] = float(data[i])
		else:
			try:
				mydict[inputs['download'][i]] = int(data[i])
			except:
				mydict[inputs['download'][i]] = str(data[i])
	with open(output_path + data[inputs['download'].index('trigger_name')] + '/' + data[inputs['download'].index('trigger_name')] + '.yaml', 'w') as f:
		yaml.dump(mydict, f)

# Determine the number of sources
def count_sources(input_path):

	n = 0
	for name in os.listdir(input_path):
		if os.path.isdir(input_path + name):
			n += 1

	return n

# Determine binning based on T90 value
def lightcurve_binning(t90):

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

# Read in available detector files, apply time/energy slices, merge and bin background
def read_detectors(detect, merge, background_time, source_time, energy_range):

	bkgd_list = []               # empty array to hold/manipulate background data

	for i in range(len(detect)):
		tte = TTE.open(detect[i]) # data
		bkgd = tte.slice_time(background_time) # background
		tte = tte.slice_time(source_time) # data time slice
		tte = tte.slice_energy(energy_range) # data energy slice
		bkgd_list.append(bkgd)
		merge.append(tte) # merge data

	return bkgd_list, merge

# Merge and bin data
def bin_data(merge, bin_rate, bin_rate2):

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
		time_sliced = tte_merge.slice_time([bins2[j], bins2[j+1]])
		time_diff = []
		for k in range(len(time_sliced.data.time)-1):
			time_diff.append(time_sliced.data.time[k+1] - time_sliced.data.time[k])
		avg_time_diff = sum(time_diff)/len(time_diff)
		if bins2[j+1] - bins2[j] < 1e-3:
			print(bins2[j+1], bins2[j], counts2[j], counts2[j] / (bins2[j+1] - bins2[j]))

	err = False
	if min(bin_sizes) > 2.5 or len(bin_sizes) == 1:
		print('Not finding event with Bayesian blocks. Not saving')
		err = True

	return rates, bins2, lc_data, nbins, err

# Fit and bin background
def fit_background(err, n_detect, b_detect, bkgd_list_nai, bkgd_list_bgo, bin_rate, t90_start, t90, nai_energy_range, bgo_energy_range):

	tot_bkgd = []                               # empty array for summing NaI/BGO background data

	# Merge and bin background
	if not err:
		if len(n_detect) > 0:
			bkgd_merge_nai = TTE.merge(bkgd_list_nai, force_unique=True) # merge NaI background data
			bkgd_phaii_nai = bkgd_merge_nai.to_phaii(bin_by_time, bin_rate, time_ref=0.0) # create time binned phaii data for NaI background
		if len(b_detect) > 0:
			bkgd_merge_bgo = TTE.merge(bkgd_list_bgo, force_unique=True) # merge BGO background data
			bkgd_phaii_bgo = bkgd_merge_bgo.to_phaii(bin_by_time, bin_rate, time_ref=0.0) # create time binned phaii data for BGO background

		# Perform background fit, interpolate, integrate over energy, and sum
		src_range = [(bkgd_phaii_nai.data.tstart[0]+0.1, t90_start-0.5), (t90+1.0, bkgd_phaii_nai.data.tstop[-1]-0.1)] # NaI background time range

		try:
			backfitter_n = BackgroundFitter.from_phaii(bkgd_phaii_nai, Polynomial, time_ranges=src_range) # initialize
			backfitter_n.fit(order=2) # NaI background fit
			bkgd_n = backfitter_n.interpolate_bins(bkgd_phaii_nai.data.tstart, bkgd_phaii_nai.data.tstop) # NaI backround interpolation
			bkgd_n = bkgd_n.integrate_energy(*nai_energy_range) # NaI backround energy integration
			tot_bkgd.append(bkgd_n)
			backfitter_b = BackgroundFitter.from_phaii(bkgd_phaii_bgo, Polynomial, time_ranges=src_range) # initialize
			backfitter_b.fit(order=2) # BGO background fit
			bkgd_b = backfitter_b.interpolate_bins(bkgd_phaii_nai.data.tstart, bkgd_phaii_nai.data.tstop) # BGO backround interpolation
			bkgd_b = bkgd_b.integrate_energy(*bgo_energy_range) # BGo backround energy integration
			tot_bkgd.append(bkgd_b)
		except:
			print('Error fitting background. Not saving')
			err = True

	return tot_bkgd, err

# Subtract background from lightcurve
def subtract_background(bb_data, lc_bkgd):

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
def plot_lightcurve(plot_path, lc_data, lc_bkgd, bins, rates, data_centers, data_heights, data_widths, t90_start, t90, name, time_range):

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
	plt.xlim(time_range[0], time_range[1])
	plt.savefig(plot_path + '%s.png' % name)
	plt.clf()

	# Plot background subtracted lightcurve
	plt.bar(data_centers, data_heights, width=data_widths, color='green', alpha=1)
	plt.axvspan(t90_start, t90, color='blue', alpha=0.5)
	plt.title(name)
	plt.xlim(time_range[0], time_range[1])
	plt.savefig(plot_path + '%s_bkg_subtracted.png' % name)
	plt.clf()

# Write lightcurve to file
def write_lightcurve(output_dir, name, bb_data, nbins, data_heights):

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

	with open(output_dir + name + '_lightcurve.dat', 'w') as f:
		f.write('# DP | Time (s) | Count rate (counts/s) \n')
		f.write('# IP LinLin\n')
		for i in range(len(start_list)):
			f.write('DP ' + str(start_list[i] - start_list[0]) + ' ' + str(rate_list[i]) + '\n')
		f.write('EN')

# Create lightcurve
def make_lightcurve(inputs, source_info, input_dir, output_dir, plot_path):

	print('Making lightcurve...')

	t_centers = []

	bin_rate, bin_rate2 = lightcurve_binning(float(source_info['t90']))

	# Find all detector files in event folder
	n_detect = glob.glob(input_dir + '/*tte_n*.fit')
	b_detect = glob.glob(input_dir + '/*tte_b*.fit')

	# Reset background data
	bkgd_phaii_nai = 0
	bkgd_phaii_bgo = 0

	merge = []                          # empty array to merge detector data
	bkgd_list_nai, merge = read_detectors(n_detect, merge, inputs['background_time_range'], inputs['source_time_range'], inputs['nai_energy_range'])
	bkgd_list_bgo, merge = read_detectors(b_detect, merge, inputs['background_time_range'], inputs['source_time_range'], inputs['bgo_energy_range'])

	rates, bins2, lc_data, nbins, err = bin_data(merge, bin_rate, bin_rate2)
	tot_bkgd, err = fit_background(err, n_detect, b_detect, bkgd_list_nai, bkgd_list_bgo, bin_rate, float(source_info['t90_start']), float(source_info['t90']), inputs['nai_energy_range'], inputs['nai_energy_range'])

	if not err:
		lc_bkgd = BackgroundRates.sum_time(tot_bkgd) # sum NaI and BGO backgrounds
		bb_data = plt.hist(bins2[:-1], bins2, weights=rates, color='green', alpha=0.7)

		data_centers, data_heights, data_widths = subtract_background(bb_data, lc_bkgd)

		plot_lightcurve(plot_path, lc_data, lc_bkgd, bins2, rates, data_centers, data_heights, data_widths, source_info['t90_start'], source_info['t90'], source_info['trigger_name'], inputs['source_time_range'])
		write_lightcurve(output_dir, source_info['trigger_name'], bb_data, nbins, data_heights)

	return err

# Create spectrum
def make_spectrum(source_info, output_dir, energy_range):

	print('Making spectrum file...')

	err = False

	if source_info['flnc_best_fitting_model'] == 'flnc_band':
		band_dict = {}
		band_dict['type'] = 'Band'
		band_dict['energy_min'] = energy_range[0]
		band_dict['energy_max'] = energy_range[1]
		band_dict['alpha'] = source_info['flnc_band_alpha']
		band_dict['beta'] = source_info['flnc_band_beta']
		band_dict['ebreak'] = source_info['flnc_band_epeak']
		with open(output_dir + source_info['trigger_name'] + '_spectrum.yaml', 'w') as outfile:
			yaml.dump(band_dict, outfile)
	elif source_info['flnc_best_fitting_model'] == 'flnc_comp':
		comptonized_dict = {}
		comptonized_dict['type'] = 'Comptonized'
		comptonized_dict['energy_min'] = energy_range[0]
		comptonized_dict['energy_max'] = energy_range[1]
		comptonized_dict['index'] = source_info['flnc_comp_index']
		comptonized_dict['epeak'] = source_info['flnc_comp_epeak']
		with open(output_dir + source_info['trigger_name'] + '_spectrum.yaml', 'w') as outfile:
			yaml.dump(comptonized_dict, outfile)
	else:
		print('For now, only using events with Band function or Comptonized best spectral fits. Not saving')
		err = True

	return err


def main():
	input_file = parse_args().yaml
	inputs = read_yaml(input_file)
	input_path, output_path, plot_path = define_paths(inputs)
	n_sources = count_sources(input_path)

	n = 0
	for name in os.listdir(input_path):
		if os.path.isdir(input_path + name):
			n += 1
			print(str(n) + '/' + str(n_sources) + ': ' + name)
			source_info = read_yaml(input_path + name + '/' + name + '.yaml')
			s_error = make_spectrum(source_info, output_path, inputs['cosi_energy_range'])
			if not s_error:
				lc_error = make_lightcurve(inputs, source_info, input_path + name, output_path, plot_path)


if __name__ == "__main__":
    main()