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

		self.altitudes = u.Quantity(altitudes)
		self.earth_zeniths = earth_zeniths

		if saa_livetime is not None:
			self.saa_livetime = u.Quantity(saa_livetime)

		if exclude is not None:
			self.exclude = np.array(exclude, dtype=bool)
		else:
			self.exclude = np.zeros(len(self.times), dtype=bool)

	@property
	def attitudes(self):
		'''
		Spacecraft attitudes, computed on demand since building one per row is expensive for long orientation files.
		'''

		return [Attitude.from_axes(x=pointing[0], z=pointing[1], frame='galactic') for pointing in self.pointings]

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

		x_pointings = SkyCoord(data[:, 2]*u.deg, data[:, 1]*u.deg, frame='galactic')
		z_pointings = SkyCoord(data[:, 4]*u.deg, data[:, 3]*u.deg, frame='galactic')
		earth_zeniths = SkyCoord(data[:, 7]*u.deg, data[:, 6]*u.deg, frame='galactic')

		pointings = list(zip(x_pointings, z_pointings))
		earth_zeniths = list(earth_zeniths)

		orientation = cls(times, pointings, altitudes, earth_zeniths)

		return orientation

	def exclude_times(self, time_range):
		'''
		Exclude time range.

		Parameters
		----------
		time_range : tuple of astropy.units.quantity.Quantity
			Range of times to exclude
		'''

		self.exclude |= (self.times >= time_range[0]) & (self.times <= time_range[1])

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

		indices = np.nonzero((orientation.times >= time_range[0]) & (orientation.times <= time_range[1]))[0]

		if len(indices) == 0 and (time_range[1] - time_range[0]) <= 1. * u.s:
			indices = np.array([np.abs(orientation.times - time_range[0]).argmin()])

		sliced_orientation = cls(orientation.times[indices], 
								 [orientation.pointings[i] for i in indices], 
								 orientation.altitudes[indices], 
								 [orientation.earth_zeniths[i] for i in indices], 
								 None, 
								 np.asarray(orientation.exclude)[indices])

		return sliced_orientation

	@classmethod
	def shift(cls, orientation, file, time):
		'''
		Shift timestamps to begin at specified time, and write to file.

		Parameters
		----------
		file : pathlib.PosixPath
			Path to .ori file
		time : astropy.units.quantity.Quantity
			New start time of orientation file

		Returns
		-------
		shifted_orientation : cosiburstpy.megalib.spacecraft_orientation.SpacecraftOrientation
			Shifted orientation
		'''

		times = orientation.times - np.min(orientation.times) + time
		keep = times.value >= 0.

		saa_livetime = orientation.saa_livetime[keep] if hasattr(orientation, 'saa_livetime') else None

		shifted_orientation = cls(times[keep], 
								  [p for p, k in zip(orientation.pointings, keep) if k], 
								  orientation.altitudes[keep], 
								  [z for z, k in zip(orientation.earth_zeniths, keep) if k], 
								  saa_livetime, 
								  np.asarray(orientation.exclude)[keep])

		shifted_orientation.write_file(file)

		return shifted_orientation

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