import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from scoords import Attitude

class SpacecraftOrientation():

	def __init__(self, times, pointings, altitudes, earth_zeniths, saa_livetime=None, exclude=None):
		'''
		Define orientation of spacecraft over time.

		Parameters
		----------
		times : list of astropy.units.quantity.Quantity
			Times
		pointings : list of 2-tuple of astropy.coordinates.sky_coordinate.SkyCoord
			Spacecraft pointings
		altitudes : list of astropy.units.quantity.Quantity
			Spacecraft altitudes
		earth_zeniths : list of astropy.coordinates.sky_coordinate.SkyCoord
			Earth zeniths
		saa_livetime : list of astropy.units.quantity.Quantity, optional
			SAA livetime
		exclude : list of bool, optional
			Whether time should be excluded
		'''

		self.times = u.Quantity(times)
		self.pointings = pointings

		self.attitudes = []
		for pointing in pointings:
			self.attitudes.append(Attitude.from_axes(x=pointing[0], z=pointing[1], frame='galactic'))

		self.altitudes = u.Quantity(altitudes)
		self.earth_zeniths = earth_zeniths

		if saa_livetime:
			self.saa_livetime = u.Quantity(saa_livetime)

		if exclude:
			self.exclude = exclude
		else:
			self.exclude = list(np.zeros(len(times), dtype=bool))

	@classmethod
	def from_file(cls, file):
		'''
		Read in orientation file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .ori file

		Returns
		-------
		orientation : cosiburstpy.spacecraft_orientation
			Spacecraft orientation
		'''

		data = np.loadtxt(file, usecols=(1, 2, 3, 4, 5, 6, 7, 8), delimiter=' ', skiprows=1, comments=("#", "EN"))

		times = data[:, 0] * u.s 
		altitudes = data[:, 5] * u.km

		pointings = []
		earth_zeniths = []

		for i in range(len(data)):

			x_pointing = SkyCoord(data[:, 2][i]*u.deg, data[:, 1][i]*u.deg, frame='galactic')
			z_pointing = SkyCoord(data[:, 4][i]*u.deg, data[:, 3][i]*u.deg, frame='galactic')
			pointings.append((x_pointing, z_pointing))

			earth_zenith = SkyCoord(data[:, 7][i]*u.deg, data[:, 6][i]*u.deg, frame='galactic')
			earth_zeniths.append(earth_zenith)

		orientation = cls(list(times), pointings, list(altitudes), earth_zeniths)

		return orientation

	def exclude_times(self, time_range):
		'''
		Exclude time range.

		Parameters
		----------
		time_range : tuple of astropy.units.quantity.Quantity
			Range of times to exclude
		'''

		for i, t in enumerate(self.times):
			if time_range[0] <= t <= time_range[1]:
				self.exclude[i] = True

	def get_orientation_at_time(self, time):
		'''
		Get the spacecraft orientation at a given time.

		Parameters
		----------
		time : astropy.units.quantity.Quantity
			Time

		Returns
		-------
		this_orientation : tuple
			Spacecraft orientation at the given time in the form (time, pointing, altitude, earth_zenith, exclude)
		'''

		if not np.min(self.times) <= time <= np.max(self.times):
			raise RuntimeError(f'Provided time ({time}) is outside the bounds of the times in the orientation file ({np.min(self.times)}, {np.max(self.times)}).')

		index = np.abs(self.times - time).argmin()

		time = self.times[index]
		pointing = self.pointings[index]
		altitude = self.altitudes[index]
		earth_zenith = self.earth_zeniths[index]
		exclude = self.exclude[index]

		this_orientation = (time, pointing, altitude, earth_zenith, exclude)

		return this_orientation

	@classmethod
	def slice(cls, orientation, time_range):
		'''
		Slice time range.

		Parameters
		----------
		orientation : cosiburstpy.megalib.spacecraft_orientation.SpacecraftOrientation
			Original orientation
		time_range : 2-tuple of astropy.units.quantity.Quantity
			Range of times

		Returns
		-------
		sliced_orientation : cosiburstpy.megalib.spacecraft_orientation.SpacecraftOrientation
			Sliced orientation
		'''

		times = []
		pointings = []
		altitudes = []
		earth_zeniths = []
		exclude = []

		for i, t in enumerate(orientation.times):

			if t >= time_range[0] and t <= time_range[1]:

				times.append(t)
				pointings.append(orientation.pointings[i])
				altitudes.append(orientation.altitudes[i])
				earth_zeniths.append(orientation.earth_zeniths[i])
				exclude.append(orientation.exclude[i])

		if len(times) == 0 and (time_range[1] - time_range[0]) <= 1. * u.s:

			closest_time = min(orientation.times, key=lambda x: abs(x - time_range[0]))
			index = np.where(orientation.times == closest_time)[0][0]

			times.append(closest_time)
			pointings.append(orientation.pointings[index])
			altitudes.append(orientation.altitudes[index])
			earth_zeniths.append(orientation.earth_zeniths[index])
			exclude.append(orientation.exclude[index])

		sliced_orientation = cls(times, pointings, altitudes, earth_zeniths, exclude)

		return sliced_orientation

	def edit_times(self, file, time):
		'''
		Update timestamps to begin at specified time, and write to file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .ori file
		time : astropy.units.quantity.Quantity
			New start time of orientation file
		'''

		start_time = np.min(self.times)
		time_add = time - start_time

		for i in range(len(self.times)):
			self.times[i] += time_add

		self.write_file(file)

	def write_file(self, file):
		'''
		Write spacecraft orientation to .ori file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .ori file
		'''

		with open(file, 'w') as f:

			f.write('# Type OrientationsGalactic\n')
			f.write('#\n')
			f.write('#\n')

			for i in range(len(self.times)):

				if hasattr(self, 'saa_livetime'):

					f.write(f'OG {self.times[i].to(u.s).value} {self.pointings[i][0].b.degree} {self.pointings[i][0].l.degree} {self.pointings[i][1].b.degree} {self.pointings[i][1].l.degree} {self.altitudes[i].to(u.km).value} {self.earth_zeniths[i].b.degree} {self.earth_zeniths[i].l.degree} {self.saa_livetime[i].to(u.s).value}\n')

				else:

					f.write(f'OG {self.times[i].to(u.s).value} {self.pointings[i][0].b.degree} {self.pointings[i][0].l.degree} {self.pointings[i][1].b.degree} {self.pointings[i][1].l.degree} {self.altitudes[i].to(u.km).value} {self.earth_zeniths[i].b.degree} {self.earth_zeniths[i].l.degree}\n')

			f.write('# EN')