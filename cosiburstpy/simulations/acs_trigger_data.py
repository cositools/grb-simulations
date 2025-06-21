import os
import csv
import h5py
from .plotting import plot_lightcurve
from .util import read_hdf5, read_csv

def plot_triggers(path, plot_path, colors=['red', 'green', 'blue', 'orange', 'purple', 'pink']):
	"""
	Plot lightcurve for each event.

	Parameters
	----------
	path : str
		Path to trigger algorithm input directory
	plot_path : str
		Path to directory to save plots
	bin_size : float, optional
		Time bin size in s
	time_range : tuple of floats, optional
		Start and end times  in s
	"""

	times = {}

	for file in os.listdir(path):

		if file.endswith('hdf5'):

			key = file[:5]

			[time_array, energy_array] = read_hdf5(path + file)[0]['trigger_data']

			mask = (energy_array >= 80.) & (energy_array <= 2000.)

			times[key] = time_array[mask]

		elif file.endswith('txt'):

			data = read_csv(path + file)

	for i in range(len(data['Event name'])):
		plot_range = [float(data['Start time (s)'][i]) - 30, float(data['Start time (s)'][i]) + 40]
		plot_lightcurve(times, save_path=plot_path + data['Event name'][i] + '.png', time_range=plot_range, event_time_range=(start_time, end_time), color=colors, title=data['Event name'][i] + ' (T90 = ' + data['GBM T90 (s)'][i] + ' s)', legend=True)