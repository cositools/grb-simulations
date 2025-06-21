import numpy as np
import matplotlib.pyplot as plt

def plot_unbinned_lightcurve(times, save_path=None, bin_size=0.05, time_range=None, event_time_range=None, color='blue', title=None, legend=False, show=False, dpi=350):
	'''
	Plot lightcurve from unbinned counts.

	Parameters
	----------
	times : list or dict of lists
		Times of events (s)
	save_path : str, optional
		Path to file to save plot
	bin_size : float, optional
		Time bin size (s)
	time_range : tuple of floats, optional
		Start and end times of plot (s)
	event_time_range : tuple of floats, optional
		Start and end times of event to shade on plot (s)
	color : str or list, optional
		Color(s) of lightcurve
	title : str, optional
		Plot title
	legend : bool, optional
		Whether to include legend on plot
	show : bool, optional
		Whether to show plot
	dpi : int, optional
		Figure resolution
	'''

	if time_range:

		nbins = round((time_range[1] - time_range[0]) / bin_size)

	if isinstance(times, list):

		if time_range:

			mask = (times >= time_range[0]) & (times <= time_range[1])
			times = times[mask]
			nbins = round((time_range[1] - time_range[0]) / bin_size)

		else:

			nbins = round((max(times) - min(times)) / bin_size)

		plot_time_range = (min(times), min(times) + (bin_size * nbins))
		bin_edges = np.linspace(plot_time_range[0], plot_time_range[1], nbins + 1)

		counts, bins = np.histogram(times, bins=bin_edges)
		plt.stairs(counts, bins, color=color, alpha=0.4)

	elif isinstance(times, dict):

		if time_range:

			nbins = round((time_range[1] - time_range[0]) / bin_size)

			for key in times:

				mask = (times[key] >= time_range[0]) & (times[key] <= time_range[1])
				times[key] = times[key][mask]

		else:

			nbins = round((max(max(t) for t in times.values()) - min(min(t) for t in times.values())) / bin_size)

		plot_time_range = (min(min(t) for t in times.values()), min(min(t) for t in times.values()) + (bin_size * nbins))
		bin_edges = np.linspace(plot_time_range[0], plot_time_range[1], nbins + 1)

		for i, key in enumerate(times):

			if isinstance(color, str):
				this_color = color
			else:
				this_color = color[i]

			counts, bins = np.histogram(times[key], bins=bin_edges)
			plt.stairs(counts, bins, color=this_color, alpha=0.4, label=key)

		if legend:
			plt.legend()

	else:

		raise RuntimeError("'event_times' must be either a list or dictionary.")

	plt.xlabel('Time (s)')
	plt.ylabel(f'Counts per {bin_size} s bin')

	if event_time_range:
		plt.axvspan(event_time_range[0], event_time_range[1], facecolor='gray', alpha=0.5)

	if title:
		plt.title(title)

	if plot_path:
		plt.savefig(save_path)

	if show:
		plt.show()

	plt.clf()
