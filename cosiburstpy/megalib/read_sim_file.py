import ROOT as root
import logging
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

	times = {'b1': [], 'b2': [], 'x1': [], 'x2': [], 'y1': [], 'y2': []}
	energies = {'b1': [], 'b2': [], 'x1': [], 'x2': [], 'y1': [], 'y2': []}

	megalib = LoadMEGAlib(mass_model)
	megalib.open_file(file)

	with SuppressOutput():

		while True:

			event = megalib.reader.GetNextEvent()

			if not event:
				break

			root.SetOwnership(event, True)

			time = float(event.GetTime().GetAsSeconds()) 

			bottom_Zplus_1 = 0.
			bottom_Zplus_2 = 0.
			bottom_Zplus_3 = 0.
			bottom_Zplus_4 = 0.
			bottom_Zplus_5 = 0.

			bottom_Zminus_1 = 0.
			bottom_Zminus_2 = 0.
			bottom_Zminus_3 = 0.
			bottom_Zminus_4 = 0.
			bottom_Zminus_5 = 0.

			x1_1 = 0.
			x1_2 = 0.
			x1_3 = 0.

			x2_1 = 0.
			x2_2 = 0.
			x2_3 = 0.

			y1_1 = 0.
			y1_2 = 0.
			y1_3 = 0.

			y2_1 = 0.
			y2_2 = 0.
			y2_3 = 0.

			for i in range(event.GetNHTs()):

				hit = event.GetHTAt(i)

				if hit.GetDetectorType() == 8:

					x = hit.GetPosition().X()
					y = hit.GetPosition().Y()
					z = hit.GetPosition().Z()

					if z < 7:

						if y > 3:

							if x < -11.2:
								bottom_Zplus_1 += hit.GetEnergy()
							elif x > -11.2 and x < -2.5:
								bottom_Zplus_2 += hit.GetEnergy()
							elif x > -2.5 and x < 4.6:
								bottom_Zplus_3 += hit.GetEnergy()
							elif x > 4.6 and x < 11.2:
								bottom_Zplus_4 += hit.GetEnergy()
							elif x > 4.6:
								bottom_Zplus_5 += hit.GetEnergy()

						elif y < 3:
               
							if x < -11.2:
								bottom_Zminus_1 += hit.GetEnergy()
							elif x > -11.2 and x < -2.5:
								bottom_Zminus_2 += hit.GetEnergy()
							elif x > -2.5 and x < 4.6:
								bottom_Zminus_3 += hit.GetEnergy()
							elif x > 4.6 and x < 11.2:
 								bottom_Zminus_4 += hit.GetEnergy()
							elif x > 4.6:
								bottom_Zminus_5 += hit.GetEnergy()

					elif z > 7:

						if x > 15:

							if y < -2.6 and y > -14:
								y1_1 += hit.GetEnergy()
							elif y > -2.6 and y < 9:
								y1_2 + hit.GetEnergy()
							elif y > 9 and y < 20.6:
								y1_3 += hit.GetEnergy()

						elif x < -10:

							if y < -2.6 and y > -14:
								y2_1 += hit.GetEnergy()
							elif y > -2.6 and y < 9:
								y2_2 += hit.GetEnergy()
							elif y > 9 and y < 20.6:
								y2_3 += hit.GetEnergy()

					elif z > -10:

						if y > 15:

							if x < -6:
								x1_1 += hit.GetEnergy()
							elif x > -6 and x < 6:
								x1_2 += hit.GetEnergy()
							elif x > 6:
								x1_3 += hit.GetEnergy()

						elif y < -10:

							if x < -6:
								x2_1 += hit.GetEnergy()
							elif x > -6 and x < 6:
								x2_2 += hit.GetEnergy()
							elif x > 6:
								x2_3 += hit.GetEnergy()

					else:

						logger.warning(f"Coordinate ({x}, {y}, {z}) not found.")

			if bottom_Zplus_1 >= 80.:
				times['b1'].append(float(time)) * u.s
				energies['b1'].append(float(bottom_Zplus_1)) * u.kev

			if bottom_Zplus_2 >= 80.:
				times['b1'].append(float(time)) * u.s
				energies['b1'].append(float(bottom_Zplus_2)) * u.kev

			if bottom_Zplus_3 >= 80.:
				times['b1'].append(float(time)) * u.s
				energies['b1'].append(float(bottom_Zplus_3)) * u.kev

			if bottom_Zplus_4 >= 80.:
				times['b1'].append(float(time)) * u.s
				energies['b1'].append(float(bottom_Zplus_4)) * u.kev

			if bottom_Zplus_5 >= 80.:
				times['b1'].append(float(time)) * u.s
				energies['b1'].append(float(bottom_Zplus_5)) * u.kev

			if bottom_Zminus_1 >= 80.:
				times['b2'].append(float(time)) * u.s
				energies['b2'].append(float(bottom_Zminus_1)) * u.kev

			if bottom_Zminus_2 >= 80.:
				times['b2'].append(float(time)) * u.s
				energies['b2'].append(float(bottom_Zminus_2)) * u.kev

			if bottom_Zminus_3 >= 80.:
				times['b2'].append(float(time)) * u.s
				energies['b2'].append(float(bottom_Zminus_3)) * u.kev

			if bottom_Zminus_4 >= 80.:
				times['b2'].append(float(time)) * u.s
				energies['b2'].append(float(bottom_Zminus_4)) * u.kev

			if bottom_Zminus_5 >= 80.:
				times['b2'].append(float(time)) * u.s
				energies['b2'].append(float(bottom_Zminus_5)) * u.kev

			if x1_1 >= 80.:
				times['x1'].append(float(time)) * u.s
				energies['x1'].append(float(x1_1)) * u.kev

			if x1_2 >= 80.:
				times['x1'].append(float(time)) * u.s
				energies['x1'].append(float(x1_2)) * u.kev

			if x1_3 >= 80.:
				times['x1'].append(float(time)) * u.s
				energies['x1'].append(float(x1_3)) * u.kev

			if x2_1 >= 80.:
				times['x2'].append(float(time)) * u.s
				energies['x2'].append(float(x2_1)) * u.kev

			if x2_2 >= 80.:
				times['x2'].append(float(time)) * u.s
				energies['x2'].append(float(x2_2)) * u.kev

			if x2_3 >= 80.:
				times['x2'].append(float(time)) * u.s
				energies['x2'].append(float(x2_3)) * u.kev

			if y1_1 >= 80.:
				times['y1'].append(float(time)) * u.s
				energies['y1'].append(float(y1_1)) * u.kev

			if y1_2 >= 80.:
				times['y1'].append(float(time)) * u.s
				energies['y1'].append(float(y1_2)) * u.kev

			if y1_3 >= 80.:
				times['y1'].append(float(time)) * u.s
				energies['y1'].append(float(y1_3)) * u.kev

			if y2_1 >= 80.:
				times['y2'].append(float(time)) * u.s
				energies['y2'].append(float(y2_1)) * u.kev

			if y2_2 >= 80.:
				times['y2'].append(float(time)) * u.s
				energies['y2'].append(float(y2_2)) * u.kev

			if y2_3 >= 80.:
				times['y2'].append(float(time)) * u.s
				energies['y2'].append(float(y2_3)) * u.kev

	acs_data = ACSData({key: list(zip(times[key], energies[key])) for key in times})

	return acs_data
	