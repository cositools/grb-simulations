from pathlib import Path
import shutil
import os
from astropy.coordinates import SkyCoord
import astropy.units as u
from scoords import Attitude, SpacecraftFrame
import numpy as np
import random
import copy

import logging
logging.basicConfig(level=logging.INFO)

from cosiburstpy import read_yaml, parse_args, write_readme, SpacecraftOrientation, Lightcurve, SourceFile, Event, Spectrum, EventTable, read_csv

input_file = parse_args([['-y', '--yaml', 'Path to input .yaml file']]).yaml
inputs = read_yaml(input_file)

input_path = Path(inputs['paths']['input'])
output_path = Path(inputs['paths']['output'])
mass_model_path = Path(inputs['paths']['mass_model'])
orientation_path = Path(inputs['paths']['orientation'])
exclusion_times_paths = [Path(inputs['paths']['exclusion_times'][0]), inputs['paths']['exclusion_times'][1]]
gbm_data_path = Path(inputs['paths']['gbm_data'])

spectral_models = inputs['events']['spectral_models']
time_between_events = [float(inputs['events']['time_separation'][0]) * u.s, float(inputs['events']['time_separation'][1]) * u.s]
energy_range = [float(inputs['events']['energy'][0]) * u.keV, float(inputs['events']['energy'][1]) * u.keV]
flux_range = [float(inputs['events']['peak_flux'][0]) / u.cm**2 / u.s, float(inputs['events']['peak_flux'][1]) / u.cm**2 / u.s]
zenith_range = [float(inputs['events']['zenith'][0]) * u.deg, float(inputs['events']['zenith'][1]) * u.deg]
azimuth_range = [float(inputs['events']['azimuth'][0]) * u.deg, float(inputs['events']['azimuth'][1]) * u.deg]
min_peak_flux = float(inputs['events']['min_peak_flux']) / u.cm**2 / u.s
start_time = float(inputs['events'].get('start_time', 0.)) * u.s

save_acs_hits = inputs['source_files']['save_acs_hits']
earth_occultation = inputs['source_files']['earth_occultation']

(output_path / 'inputs').mkdir(parents=True, exist_ok=True)
(output_path / 'source_files').mkdir(parents=True, exist_ok=True)
shutil.copy(os.path.abspath(__file__), output_path / 'inputs')
shutil.copy(input_file, output_path / 'inputs')

write_readme(output_path / 'README.md', inputs_path=output_path / 'inputs', input_parameters=inputs)

events = {}

for spectrum_file in (input_path / 'spectra').iterdir():

	if spectrum_file.suffix == '.yaml':

		spectrum = Spectrum.from_file(input_path / 'spectra' / spectrum_file)

		if spectrum.model in spectral_models:

			name = spectrum_file.name[:11]

			gbm_data = read_yaml(gbm_data_path / name / f'{name}.yaml')
			peak_flux_model = gbm_data['pflx_best_fitting_model']

			if not peak_flux_model == 'nan':
				peak_flux = gbm_data[peak_flux_model + '_phtfluxb'] / u.cm**2 / u.s

				if peak_flux > min_peak_flux: 

					lightcurve_file = None

					for item in (input_path / 'lightcurves').iterdir():

						if name in str(item) and item.suffix == '.dat':
							lightcurve_file = item

					if lightcurve_file:

						lightcurve = Lightcurve.from_file(lightcurve_file)
						events[name] = (spectrum, lightcurve, lightcurve_file)

events = list(events.items())
random.shuffle(events)
events = dict(events)

exclusion_times = []

for file in exclusion_times_paths:

	data = read_csv(file, delimiter=',')
	data = {key: [float(item) for item in entry] for key, entry in data.items()}

	time_ranges = list(zip(data['start'] * u.s, data['stop'] * u.s))
	exclusion_times.extend(time_ranges)

event_table_dict = {'Event': [], 'Start time (s)': [], 'Peak start time (s)': [], 'Lightcurve duration (s)': [], 'GBM T90 (s)': [], 
					'Peak photon flux (ph/cm^2/s) for 10-10000 keV': [], 'Peak energy flux (erg/cm^2/s) for 10-10000 keV': [], 
					'Peak photon flux (ph/cm^2/s) for 80-2000 keV': [], 'Peak photon flux (ph/cm^2/s) for 50-300 keV': [], 
					'Peak photon flux (ph/cm^2/s) for 10-1000 keV': [], 'Zenith (deg)': [], 'Azimuth (deg)': [], 
					'Galactic latitude (deg)': [], 'Galactic longitude (deg)': [], 'Spectral model': [], 'Spectral parameters': []}

orientation = SpacecraftOrientation.from_file(orientation_path)

for time_range in exclusion_times:
	orientation.exclude_times(time_range)

time = orientation.times.min() + start_time

for name in events:

	source_file_path = output_path / 'source_files' / f'{name}.source'

	spectrum = events[name][0]
	lightcurve = events[name][1]
	
	duration = lightcurve.bin_edges[-1] - lightcurve.bin_edges[0]

	time_add = np.random.uniform(time_between_events[0].to(u.s).value, time_between_events[1].to(u.s).value) * u.s

	exclude = True
	while exclude == True:

		exclude_start = orientation.get_orientation_at_time(time + time_add)[4]
		exclude_end = orientation.get_orientation_at_time(time + time_add + duration)[4]

		if any((exclude_start, exclude_end)):
			time_add += np.random.uniform(time_between_events[0].to(u.s).value, time_between_events[1].to(u.s).value) * u.s

		else:
			exclude = False

	event_time = time + time_add
	time = event_time + duration

	# MEGAlib runs for the lightcurve duration + 10 s (see SourceFile), so only this window of the orientation is needed
	orientation_window = SpacecraftOrientation.slice(orientation, (event_time - 1. * u.s, event_time + duration + 11. * u.s))
	this_orientation = SpacecraftOrientation.shift(orientation_window, output_path / 'source_files' / f'{name}.ori', 
												orientation_window.times.min() - round(event_time.to(u.s).value) * u.s)

	lightcurve.edit_times(output_path / 'source_files' / f'{name}.dat', 0. * u.s)

	event_orientation = SpacecraftOrientation.slice(this_orientation, (lightcurve.bin_edges[0], lightcurve.bin_edges[-1]))

	attitude = Attitude.from_axes(x=event_orientation.pointings[0][0], z=event_orientation.pointings[0][1], frame='galactic')

	peak_flux = np.exp(np.random.uniform(np.log(flux_range[0].value), np.log(flux_range[1].value))) * flux_range[0].unit

	event = Event(name, lightcurve, spectrum, event_orientation)
	event.set_peak_flux(peak_flux, energy_range)

	occulted = True
	while occulted == True:

		sin_dec_range = [np.sin(np.radians((90.*u.deg - zenith_range[1]).to(u.deg).value)), np.sin(np.radians((90.*u.deg - zenith_range[0]).to(u.deg).value))]
		sin_dec = np.random.uniform(sin_dec_range[0], sin_dec_range[1])
		dec = np.degrees(np.arcsin(sin_dec)) * u.deg

		zenith = 90. * u.deg - dec
		azimuth = np.random.uniform(azimuth_range[0].to(u.deg).value, azimuth_range[1].to(u.deg).value) * u.deg

		position = SkyCoord(azimuth, 90.*u.deg - zenith, representation_type='spherical', frame=SpacecraftFrame(attitude=attitude))
		event.position = position

		occulted = event.occulted

	source_file = SourceFile(event, mass_model_path, events[name][2].name, Path(f'{name}.ori'))
	source_file.write_file(source_file_path, save_acs_hits, earth_occultation=earth_occultation)

	gbm_data = read_yaml(gbm_data_path / name / f'{name}.yaml')

	event_table_dict['Event'].append(name)
	event_table_dict['Start time (s)'].append(f'{event_time.to(u.s).value:.2f}')
	event_table_dict['Peak start time (s)'].append(f'{lightcurve.peak[0].to(u.s).value:.2f}')
	event_table_dict['Lightcurve duration (s)'].append(f'{duration.to(u.s).value:.2f}')
	event_table_dict['GBM T90 (s)'].append(gbm_data['t90'])
	event_table_dict['Peak photon flux (ph/cm^2/s) for 10-10000 keV'].append(f'{event.peak_photon_flux((10.*u.keV, 10000.*u.keV)).to(u.cm**-2*u.s**-1).value:.3f}')
	event_table_dict['Peak energy flux (erg/cm^2/s) for 10-10000 keV'].append(f'{event.peak_energy_flux((10.*u.keV, 10000.*u.keV)).to(u.erg/u.cm**2/u.s).value:.3e}')
	event_table_dict['Peak photon flux (ph/cm^2/s) for 80-2000 keV'].append(f'{event.peak_photon_flux((80.*u.keV, 2000.*u.keV)).to(u.cm**-2*u.s**-1).value:.3f}')
	event_table_dict['Peak photon flux (ph/cm^2/s) for 50-300 keV'].append(f'{event.peak_photon_flux((50.*u.keV, 300.*u.keV)).to(u.cm**-2*u.s**-1).value:.4f}')
	event_table_dict['Peak photon flux (ph/cm^2/s) for 10-1000 keV'].append(f'{event.peak_photon_flux((10.*u.keV, 1000.*u.keV)).to(u.cm**-2*u.s**-1).value:.4f}')
	event_table_dict['Zenith (deg)'].append(f'{90. - event.position.lat.to(u.deg).value:.2f}')
	event_table_dict['Azimuth (deg)'].append(f'{event.position.lon.to(u.deg).value:.2f}')
	event_table_dict['Galactic latitude (deg)'].append(f'{event.position.transform_to("galactic").b.to(u.deg).value:.2f}')
	event_table_dict['Galactic longitude (deg)'].append(f'{event.position.transform_to("galactic").l.to(u.deg).value:.2f}')
	event_table_dict['Spectral model'].append(event.spectrum.model)
	event_table_dict['Spectral parameters'].append(event.spectrum.parameter_list)

keys = list(event_table_dict.keys())
sort_index = keys.index('Start time (s)')

rows = list(zip(*[event_table_dict[k] for k in keys]))
sorted_rows = sorted(rows, key=lambda row: float(row[sort_index]))

event_table_dict = {key: [row[i] for row in sorted_rows] for i, key in enumerate(keys)}

event_table = EventTable(event_table_dict)
event_table.write_file(output_path / 'events.csv')

