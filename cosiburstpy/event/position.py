import numpy as np
from astropy.coordinates import Angle
import astropy.units as u

def check_if_earth_occulted(position, earth_zenith, altitude):
	'''
	Determine if position is occulted by the Earth.

	Parameters
	----------
	position : astropy.coordinates.sky_coordinate.SkyCoord
		Position
	earth_zenith : astropy.coordinates.sky_coordinate.SkyCoord
		Earth zenith
	altitude : astropy.units.quantity.Quantity
		Altitude of spacecraft

	Returns
	-------
	occulted : bool
		Whether source is Earth occulted
	'''

	r_earth = 6378. * u.km

	source_angle = position.separation(earth_zenith)
	max_angle = Angle(np.pi * u.rad - np.arcsin(r_earth / (r_earth + altitude)))

	if source_angle.wrap_at(180*u.deg).deg >= max_angle.wrap_at(180*u.deg).deg:
		occulted = True
	else:
		occulted = False

	return occulted
