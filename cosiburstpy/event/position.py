import numpy as np
from astropy.coordinates import SkyCoord, Angle
import astropy.units as u
from scoords import Attitude, SpacecraftFrame

def galactic_from_local_position(position, pointing):
	'''
	Calculate position in galactic coordinates from spacecraft pointing, zenith, and azimuth.

	Parameters
	----------
	position : astropy.coordinates.sky_coordinate.SkyCoord
		Position in spaecraft coordinates
	pointing : tuple of astropy.coordinates.sky_coordinate.SkyCoord
		Spacecraft pointing in the form (x-axis pointing, z-axis pointing)

	Returns
	-------
	position_galactic : astropy.coordinates.sky_coordinate.SkyCoord
		Position in galactic coordinates
	'''

	attitude = Attitude.from_axes(x=pointing[0], z=pointing[1], frame='galactic')
	position_galactic = SkyCoord(position.lon, position.lat, representation_type='spherical', frame=SpacecraftFrame(attitude=attitude)).transform_to('galactic')

	return position_galactic

def local_from_galactic_position(position, pointing):
	'''
	Calculate position in local coordinates from spacecraft pointing, longitude, and latitude.

	Parameters
	----------
	position : astropy.coordinates.sky_coordinate.SkyCoord
		Position in galactic coordinates
	pointing : tuple of astropy.coordinates.sky_coordinate.SkyCoord
		Spacecraft pointing in the form (x-axis pointing, z-axis pointing)

	Returns
	-------
	position_local : astropy.coordinates.sky_coordinate.SkyCoord
		Position in spaecraft coordinates
	'''

	attitude = Attitude.from_axes(x=pointing[0], z=pointing[1], frame='galactic')
	position_local = SkyCoord(position.l, position.b, representation_type='spherical', frame='galactic').transform_to(SpacecraftFrame(attitude=attitude))

	return position_local

def check_if_earth_occulted(source_position, earth_zenith, altitude):
	'''
	Determine if source is occulted by the Earth for a particular altitude and Earth zenith.

	Parameters
	----------
	source_position : astropy.coordinates.sky_coordinate.SkyCoord
		Position of source
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

	source_angle = source_position.separation(earth_zenith)
	max_angle = Angle(np.pi * u.rad - np.arcsin(r_earth / (r_earth + altitude)))

	if source_angle.wrap_at(180*u.deg).deg >= max_angle.wrap_at(180*u.deg).deg:
		occulted = True
	else:
		occulted = False

	return occulted







