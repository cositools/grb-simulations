import argparse
import yaml
import os
import numpy as np

# Initialize parser & read arguments from command line
def parse_args():

	parser = argparse.ArgumentParser()
	parser.add_argument("-y", "--yaml", help = "Path to input .yaml file", default='examples/example.yaml')
	args = parser.parse_args()

	return args

# Read yaml file
def read_yaml(file):

	with open(file, "r") as myfile:
		inputs = yaml.safe_load(myfile)

	return inputs

# Define paths to input files, output source files, and current working directory
def define_paths(inputs):

	if 'event_subtype' in inputs:
		input_path = os.path.abspath(inputs['input_path'] + inputs['event_type'] + '/' + inputs['event_subtype']) + '/'
		output_path = os.path.abspath(inputs['output_path'] + inputs['event_type'] + '/' + inputs['event_subtype']) + '/'
	else:
		input_path = os.path.abspath(inputs['input_path'] + inputs['event_type']) + '/'
		output_path = os.path.abspath(inputs['output_path'] + inputs['event_type']) + '/'
	if not os.path.isdir(output_path):
		os.makedirs(output_path)

	return input_path, output_path

# Create list of event names and directories with spectrum and lightcurve file names
def make_event_dict(input_path, mix_or_match, spectrum_type):

	key_list = []
	lightcurve_dict = {}
	spectrum_dict = {}
	for file in os.listdir(input_path):
		filename = os.fsdecode(file)
		if filename.endswith('_lightcurve.dat'):
			key_list.append(filename[:-15])
			lightcurve_dict[filename[:-15]] = filename

	if mix_or_match == 'match':
		if spectrum_type == 'yaml':
			for item in key_list:
				filename = item + '_spectrum.yaml'
				spectrum_dict[item] = filename
		elif spectrum_type == 'dat':
			for item in key_list:
				filename = item + '_spectrum.dat'
				spectrum_dict[item] = filename
		else:
			raise RuntimeError("'spectrum_type' must be either 'yaml' or 'dat'.")
	elif mix_or_match == 'mix':
		spectrum_list = []
		for file in os.listdir(input_path):
			filename = os.fsdecode(file)
			if spectrum_type == 'yaml':
				if filename.endswith('_spectrum.yaml'):
					spectrum_list.append(filename)
			elif spectrum_type == 'dat':
				if filename.endswith('_spectrum.dat'):
					spectrum_list.append(filename)
			else:
				raise RuntimeError("'spectrum_type' must be either 'yaml' or 'dat'.")
		for item in key_list:
			index = np.random.randint(len(spectrum_list))
			filename = spectrum_list[index]
			spectrum_dict[item] = filename
	else:
		raise RuntimeError("'mix_or_match' must be either 'mix' or 'match'.")

	return key_list, lightcurve_dict, spectrum_dict

# Defines zenith, azimuth, and flux of source
def define_angles_flux(z_min, z_max, a_min, a_max, flux_min, flux_max):

	zenith = np.random.uniform(float(z_min), float(z_max))
	azimuth = np.random.uniform(float(a_min), float(a_max))
	flux = np.random.uniform(float(flux_min), float(flux_max))

	return zenith, azimuth, flux

# Determines text to define spectrum in source file
def define_spectrum(spectrum, name):

	if spectrum['type'] == 'File':
		spectrum_text = name + '.Spectrum        File ' + spectrum['filename']
	elif spectrum['type'] == 'Mono':
		spectrum_text = name + '.Spectrum        Mono ' + str(spectrum['energy'])
	elif spectrum['type'] == 'Band':
		spectrum_text = name + '.Spectrum        Band ' + str(spectrum['energy_min']) + ' ' + str(spectrum['energy_max']) + ' ' + str(spectrum['alpha']) + ' ' + str(spectrum['beta']) + ' ' + str(spectrum['ebreak'])
	elif spectrum['type'] == 'Comptonized':
		spectrum_text = name + '.Spectrum        Comptonized ' + str(spectrum['energy_min']) + ' ' + str(spectrum['energy_max']) + ' ' + str(spectrum['index']) + ' ' + str(spectrum['epeak'])
	elif spectrum['type'] == 'PowerLaw':
		spectrum_text = name + '.Spectrum        PowerLaw ' + str(spectrum['energy_min']) + ' ' + str(spectrum['energy_max']) + ' ' + str(spectrum['index'])
	elif spectrum['type'] == 'BrokenPowerLaw':
		spectrum_text = name + '.Spectrum        BrokenPowerLaw ' + str(spectrum['energy_min']) + ' ' + str(spectrum['energy_max']) + ' ' + str(spectrum['ebreak']) + ' ' + str(spectrum['index_lo']) + ' ' + str(spectrum['index_hi'])
	else:
		raise RuntimeError("Spectral type not supported. 'type' in spectrum .yaml file must be 'Mono', 'Band', 'Comptonized', 'PowerLaw', or 'BrokenPowerLaw'.")

	return spectrum_text

# Reads spectra and outputs text for lightcurve and spectrum in source file
def lightcurve_spectrum_text(name, lightcurve, event_spectrum, path):

	lightcurve_text = name + '.Lightcurve      File false ' + path + lightcurve

	if event_spectrum.endswith('_spectrum.yaml'):
		spectrum = read_yaml(path + event_spectrum)
	elif event_spectrum.endswith('_spectrum.dat'):
		spectrum = {}
		spectrum['type'] = 'File'
		spectrum['filename'] = path + event_spectrum
	else:
		raise RuntimeError("Problem with spectrum file name.")

	spectrum_text = define_spectrum(spectrum, name)

	return spectrum_text, lightcurve_text

# Creates source file
def make_source_file(event, filename, geometry, shield_counts, z, a, flux, spectrum_text, lightcurve_text):

	with open(filename, 'w') as f:
		f.write('# Global Parameters\n')
		f.write('Version                     1\n')
		f.write('Geometry                    ' + geometry + '\n')
		f.write('\n# Physics list\n')  
		f.write('PhysicsListEM               LivermorePol\n')
		f.write('\n# Output formats\n')
		f.write('StoreSimulationInfo         init-only\n')
		if shield_counts == 'y':
			f.write('\n# Store shield counts\n')
			f.write('PreTriggerMode              EveryEventWithHits\n')
		f.write('\n# Run and source parameters\n')
		f.write('Run                         GRBSim\n')
		f.write('GRBSim.FileName             ' + event + '\n')
		f.write('GRBSim.Time                 5000.0\n')
		f.write('GRBSim.Source               ' + event + '\n')
		f.write(event + '.ParticleType    1\n')
		f.write(event + '.Beam            FarFieldPointSource ' + str(z) + ' ' + str(a) + '\n')
		f.write('\n# Spectrum \n')
		f.write(spectrum_text + '\n')
		f.write('\n# Average photon flux in photon/cm2/s\n')
		f.write(event + '.Flux            ' + str(flux) + '\n')
		f.write('\n# Lightcurve\n')
		f.write(lightcurve_text)


def main():
	input_file = parse_args().yaml
	inputs = read_yaml(input_file)
	input_path, output_path = define_paths(inputs)
	events, lightcurves, spectra = make_event_dict(input_path, inputs['mix_or_match'], inputs['spectrum_type'])
	for event in events:
		source_filename = output_path + event + '.source'
		zenith_angle, azimuthal_angle, flux = define_angles_flux(inputs['incidence_angle_min'], inputs['incidence_angle_max'], 0, 360, inputs['flux_min'], inputs['flux_max'])
		spectrum_text, lightcurve_text = lightcurve_spectrum_text(event, lightcurves[event], spectra[event], input_path)
		make_source_file(event, source_filename, inputs['geometry_path'], inputs['shield_counts'], zenith_angle, azimuthal_angle, flux, spectrum_text, lightcurve_text)

if __name__ == "__main__":
    main()
