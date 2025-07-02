import numpy as np
import astropy.units as u
import logging

logger = logging.getLogger(__name__)

class Lightcurve():

	def __init__(self, times, amplitudes):
		'''
		Define lightcurve of event.

		Parameters
		----------
		times : numpy.ndarray of astropy.units.quantity.Quantity, shape (N,)
			Times
		amplitudes : numpy.ndarray of astropy.units.quantity.Quantity, shape (N,)
			Lightcurve shape in units equivalent to counts/s
		'''

		self.bin_edges = times
		self.shape = amplitudes[:-1]

		self.peak = (self.bin_edges[np.where(self.shape == np.max(self.shape))[0][0]], self.bin_edges[np.where(self.shape == np.max(self.shape))[0][-1]+1])
		self.peak_duration = self.peak[1] - self.peak[0]

	@classmethod
	def from_file(cls, file):
		'''
		Read in lightcurve file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .dat file

		Returns
		-------
		lightcurve : cosiburstpy.lightcurve.Lightcurve
			Lightcurve
		'''

		data = np.loadtxt(file, usecols=(1, 2), delimiter=' ', skiprows=1, comments=("#", "EN"))

		times = data[:, 0] * u.s 
		amplitudes = data[:, 1] / u.s

		lightcurve = cls(times, amplitudes)

		return lightcurve

	def edit_times(self, file, time):
		'''
		Update timestamps to begin at specified time, and write to file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .dat file
		time : astropy.units.quantity.Quantity
			New start time of lightcurve
		'''

		start_time = np.min(self.bin_edges)
		time_add = time - start_time

		for i in range(len(self.bin_edges)):
			self.bin_edges[i] += time_add

		self.peak = (self.peak[0] + time_add, self.peak[1] + time_add)
		self.peak_duration = self.peak[1] - self.peak[0]

		self.write_file(file)

	def write_file(self, file):
		'''
		Write lightcurve to .dat file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .dat file
		'''

		logger.info(f"Writing file: {file}")

		(file.parent).mkdir(parents=True, exist_ok=True)

		with open(file, 'w') as f:

			f.write('# DP | Time (s) | Count rate (counts/s) \n')
			f.write('# IP LinLin\n')

			for i in range(len(self.shape)):
				f.write(f'DP {self.bin_edges[i].to(u.s).value} {self.shape[i].to(u.s**-1).value}\n')
				
			f.write(f'DP {self.bin_edges[-1].to(u.s).value} {self.shape[-1].to(u.s**-1).value}\n')
			f.write('EN')

	@property
	def peak_ratio(self):
		'''
		Ratio of peak amplitude to average amplitude of lightcurve.

		Returns
		-------
		ratio : astropy.units.quantity.Quantity
			Peak ratio
		'''

		bin_sizes = np.diff(self.bin_edges)
		ratio = np.max(self.shape) / np.dot(bin_sizes, self.shape) * np.sum(bin_sizes)

		return ratio
