import ROOT as root
import logging
import astropy.units as u
from cosiburstpy.utility.utility import SuppressOutput
from cosiburstpy.simulations.acs_data import ACSData
from .load_megalib import LoadMEGAlib

logger = logging.getLogger(__name__)

def read_sim_file(file, mass_model):
	'''
	Extract ACS hits from .sim or .sim.gz file. Modified from Savitri's code to create ACS data .csv files with the Data Challenge 3 mass model.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to ACS data .sim or .sim.gz file
	mass_model : pathlib.PosixPath
		Path to Data Challenge 3 analysis mass model

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

					x = hit.GetPosition().X()
					y = hit.GetPosition().Y()
					z = hit.GetPosition().Z()

					if x < -11.2 and y < 3 and z < 7:
						z0_0 += hit.GetEnergy()
					elif -11.2 < x < -2.5 and y < 3 and z < 7:
						z0_1 += hit.GetEnergy()
					elif -2.5 < x < 4.6 and y < 3 and z < 7:
						z0_2 += hit.GetEnergy()
					elif 4.6 < x < 11.2 and y < 3 and z < 7:
 						z0_3 += hit.GetEnergy()
					elif x > 11.2 and y < 3 and z < 7:
						z0_4 += hit.GetEnergy()

					elif x < -11.2 and y > 3 and z < 7:
						z1_0 += hit.GetEnergy()
					elif -11.2 < x < -2.5 and y > 3 and z < 7:
						z1_1 += hit.GetEnergy()
					elif -2.5 < x < 4.6 and y > 3 and z < 7:
						z1_2 += hit.GetEnergy()
					elif 4.6 < x < 11.2 and y > 3 and z < 7:
						z1_3 += hit.GetEnergy()
					elif x > 11.2 and y > 3 and z < 7:
						z1_4 += hit.GetEnergy()

					elif z > 7 and -14 < y < -2.6 and x > 15:
						y0_0 += hit.GetEnergy()
					elif z > 7 and -2.6 < y < 9 and x > 15:
						y0_1 + hit.GetEnergy()
					elif z > 7 and 9 < y < 20.6 and x > 15:
						y0_2 += hit.GetEnergy()

					elif z > 7 and -14 < y < -2.6 and x < -10:
						y1_0 += hit.GetEnergy()
					elif z > 7 and -2.6 < y < 9 and x < -10:
						y1_1 += hit.GetEnergy()
					elif z > 7 and 9 < y < 20.6 and x < -10:
						y1_2 += hit.GetEnergy()

					elif z > -10 and x < -6 and y < -10:
						x0_0 += hit.GetEnergy()
					elif z > -10 and -6 < x < 6 and y < -10:
						x0_1 += hit.GetEnergy()
					elif z > -10 and x > 6 and y < -10:
						x0_2 += hit.GetEnergy()

					elif z > -10 and x < -6 and y > 15:
						x1_0 += hit.GetEnergy()
					elif z > -10 and -6 < x < 6 and y > 15:
						x1_1 += hit.GetEnergy()
					elif z > -10 and x > 6 and y > 15:
						x1_2 += hit.GetEnergy()

					else:

						logger.warning(f"Coordinate ({x}, {y}, {z}) not found.")

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
	