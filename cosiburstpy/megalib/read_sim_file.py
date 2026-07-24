import ROOT as root
import logging
import astropy.units as u
from cosiburstpy.utility.utility import SuppressOutput
from cosiburstpy.simulations.acs_data import ACSData
from .load_megalib import LoadMEGAlib

logger = logging.getLogger(__name__)

def read_sim_file(file, mass_model):
	'''
	Extract ACS hits from .sim or .sim.gz file. Modified from Nicolò's code in cosipy.nonimaging.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to ACS data .sim or .sim.gz file
	mass_model : pathlib.PosixPath
		Path to mass model with labeled ACS crystals

	Returns
	-------
	acs_data : cosiburstpy.simulations.acs_data.ACSData
		ACS data
	'''

	logger.info(f"Reading file: {file}")

	times = {'z0': [], 'z1': [], 'x0': [], 'x1': [], 'y0': [], 'y1': []}
	energies = {'z0': [], 'z1': [], 'x0': [], 'x1': [], 'y0': [], 'y1': []}

	megalib = LoadMEGAlib(mass_model)
	megalib.open_file(file)
	geometry = megalib.geometry

	with SuppressOutput():

		while True:

			event = megalib.reader.GetNextEvent()

			if not event:
				break

			root.SetOwnership(event, True)

			time = float(event.GetTime().GetAsSeconds()) 

			z1_0 = 0.
			z1_1 = 0.
			z1_2 = 0.
			z1_3 = 0.
			z1_4 = 0.

			z0_0 = 0.
			z0_1 = 0.
			z0_2 = 0.
			z0_3 = 0.
			z0_4 = 0.

			x1_0 = 0.
			x1_1 = 0.
			x1_2 = 0.

			x0_0 = 0.
			x0_1 = 0.
			x0_2 = 0.

			y0_0 = 0.
			y0_1 = 0.
			y0_2 = 0.

			y1_0 = 0.
			y1_1 = 0.
			y1_2 = 0.

			for i in range(event.GetNHTs()):

				hit = event.GetHTAt(i)

				if hit.GetDetectorType() == 8:

					position = hit.GetPosition()

					detector_object = geometry.GetDetector(position)

					if not detector_object:
						raise RuntimeError(f"Coordinate ({position.X()}, {position.Y()}, {position.Z()}) not found.")
					else:	
						detector = detector_object.GetName()

					if detector.GetString() == 'ACS_Z0_0':
						z0_0 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Z0_1':
						z0_1 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Z0_2':
						z0_2 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Z0_3':
						z0_3 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Z0_4':
						z0_4 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Z1_0':
						z1_0 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Z1_1':
						z1_1 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Z1_2':
						z1_2 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Z1_3':
						z1_3 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Z1_4':
						z1_4 += hit.GetEnergy()

					elif detector.GetString() == 'ACS_Y0_0':
						y0_0 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Y0_1':
						y0_1 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Y0_2':
						y0_2 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Y1_0':
						y1_0 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Y1_1':
						y1_1 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_Y1_2':
						y1_2 += hit.GetEnergy()

					elif detector.GetString() == 'ACS_X0_0':
						x0_0 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_X0_1':
						x0_1 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_X0_2':
						x0_2 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_X1_0':
						x1_0 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_X1_1':
						x1_1 += hit.GetEnergy()
					elif detector.GetString() == 'ACS_X1_2':
						x1_2 += hit.GetEnergy()

					else:

						logger.warning(f"Coordinate ({position.X()}, {position.Y()}, {position.Z()}) not found.")

			if z0_0 >= 80.:
				times['z0'].append(float(time) * u.s)
				energies['z0'].append(float(z0_0) * u.keV)

			if z0_1 >= 80.:
				times['z0'].append(float(time) * u.s)
				energies['z0'].append(float(z0_1) * u.keV)

			if z0_2 >= 80.:
				times['z0'].append(float(time) * u.s)
				energies['z0'].append(float(z0_2) * u.keV)

			if z0_3 >= 80.:
				times['z0'].append(float(time) * u.s)
				energies['z0'].append(float(z0_3) * u.keV)

			if z0_4 >= 80.:
				times['y0'].append(float(time) * u.s)
				energies['y0'].append(float(z0_4) * u.keV)

			if z1_0 >= 80.:
				times['y1'].append(float(time) * u.s)
				energies['y1'].append(float(z1_0) * u.keV)

			if z1_1 >= 80.:
				times['z1'].append(float(time) * u.s)
				energies['z1'].append(float(z1_1) * u.keV)

			if z1_2 >= 80.:
				times['z1'].append(float(time) * u.s)
				energies['z1'].append(float(z1_2) * u.keV)

			if z1_3 >= 80.:
				times['z1'].append(float(time) * u.s)
				energies['z1'].append(float(z1_3) * u.keV)

			if z1_4 >= 80.:
				times['z1'].append(float(time) * u.s)
				energies['z1'].append(float(z1_4) * u.keV)

			if x0_0 >= 80.:
				times['x0'].append(float(time) * u.s)
				energies['x0'].append(float(x0_0) * u.keV)

			if x0_1 >= 80.:
				times['x0'].append(float(time) * u.s)
				energies['x0'].append(float(x0_1) * u.keV)

			if x0_2 >= 80.:
				times['x0'].append(float(time) * u.s)
				energies['x0'].append(float(x0_2) * u.keV)

			if x1_0 >= 80.:
				times['x1'].append(float(time) * u.s)
				energies['x1'].append(float(x1_0) * u.keV)

			if x1_1 >= 80.:
				times['x1'].append(float(time) * u.s)
				energies['x1'].append(float(x1_1) * u.keV)

			if x1_2 >= 80.:
				times['x1'].append(float(time) * u.s)
				energies['x1'].append(float(x1_2) * u.keV)

			if y0_0 >= 80.:
				times['y0'].append(float(time) * u.s)
				energies['y0'].append(float(y0_0) * u.keV)

			if y0_1 >= 80.:
				times['y0'].append(float(time) * u.s)
				energies['y0'].append(float(y0_1) * u.keV)

			if y0_2 >= 80.:
				times['y0'].append(float(time) * u.s)
				energies['y0'].append(float(y0_2) * u.keV)

			if y1_0 >= 80.:
				times['y1'].append(float(time) * u.s)
				energies['y1'].append(float(y1_0) * u.keV)

			if y1_1 >= 80.:
				times['y1'].append(float(time) * u.s)
				energies['y1'].append(float(y1_1) * u.keV)

			if y1_2 >= 80.:
				times['y1'].append(float(time) * u.s)
				energies['y1'].append(float(y1_2) * u.keV)

	acs_data = ACSData({key: list(zip(times[key], energies[key])) for key in times})

	return acs_data
	