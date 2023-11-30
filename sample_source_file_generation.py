import argparse
import yaml
import os
import numpy as np
from scipy import integrate

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

# Define paths to input files and output source files
def define_paths(inputs):

	input_path = os.path.abspath(inputs['input_path']) + '/'
	output_path = os.path.abspath(inputs['output_path']) + '/'
	if not os.path.isdir(output_path):
		os.makedirs(output_path)

	return input_path, output_path

# Create README file
def write_readme(inputs, input_path, output_path):

	with open(output_path + 'README.md', 'w') as f:
		f.write('These source files were created using the input files in ' + input_path + '. ')
		
		if inputs['coordinate_system'] == 'local':
			if 'zenith' in inputs:
				if type(inputs['zenith']) == int:
					f.write('All .source files have a zenith angle of ' + str(inputs['zenith']) + ' degrees. ')
				else:
					f.write('The zenith angles of each source were chosen from the following list (in degrees): ' + str(inputs['zenith']) + '. ')
			elif 'zenith_min' in inputs and 'zenith_max' in inputs:
				f.write('The zenith angles of each source were chosen randomly between ' + str(inputs['zenith_min']) + ' and ' + str(inputs['zenith_max']) + ' degrees. ')
			else:
				raise RuntimeError("Must specify zenith angle(s) in input .yaml file.")
			if 'azimuth' in inputs:
				if type(inputs['azimuth']) == int:
					f.write('All .source files have an azimuthal angle of ' + str(inputs['azimuth']) + ' degrees. ')
				else:
					f.write('The azimuthal angles of each source were chosen from the following list (in degrees): ' + str(inputs['azimuth']) + '. ')
			elif 'azimuth_min' in inputs and 'azimuth_max' in inputs:
				f.write('The azimuthal angles of each source were chosen randomly between ' + str(inputs['azimuth_min']) + ' and ' + str(inputs['azimuth_max']) + ' degrees. ')
			else:
				raise RuntimeError("Must specify azimuthal angle(s) in input .yaml file.")
		
		else:
			raise RuntimeError("Only detector coordinates are supported for now. 'coordinate_system' in input .yaml file must be 'local'.")

		if 'ph_flux' in inputs:
			if type(inputs['ph_flux']) == int:
				f.write('All .source files have a flux of ' + str(inputs['ph_flux']) + ' ph/cm<sup>2</sup>/s. ')
			else:
				f.write('The fluxes of each source were chosen from the following list (in ph/cm<sup>2</sup>/s): ' + str(inputs['ph_flux']) + '. ')
		elif 'ph_flux_min' in inputs and 'ph_flux_max' in inputs:
			f.write('The fluxes of each source were chosen randomly between ' + str(inputs['ph_flux_min']) + ' and ' + str(inputs['ph_flux_max']) + ' ph/cm<sup>2</sup>/s. ')
		elif 'e_flux' in inputs:
			if type(inputs['flux']) == int:
				f.write('All .source files have a flux of ' + str(inputs['e_flux']) + ' erg/cm<sup>2</sup>/s. ')
			else:
				f.write('The fluxes of each source were chosen from the following list (in erg/cm<sup>2</sup>/s): ' + str(inputs['e_flux']) + '. ')
		elif 'e_flux_min' in inputs and 'e_flux_max' in inputs:
			f.write('The fluxes of each source were chosen randomly between ' + str(inputs['e_flux_min']) + ' and ' + str(inputs['e_flux_max']) + ' erg/cm<sup>2</sup>/s. ')
		else:
			raise RuntimeError("Must specify flux(es) in input .yaml file.")

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

# Band function
def band_function(e, alpha, beta, ebreak):

	if e <= (alpha - beta) * ebreak:
		return (e / 100)**alpha * np.exp(-e / ebreak)
	else:
		return (e / 100)**beta * np.exp(beta - alpha) * ((alpha - beta) * ebreak / 100)**(alpha - beta)

# Comptonized function
def comp_function(e, index, epeak):

	return e**index * np.exp(-(index + 2) * e / epeak)

# Power law function
def powerlaw_function(e, index):

	return e**(-index)

# Broken power law function
def brokenpowerlaw_function(e, ebreak, index_lo, index_hi, e_max):

	if e <= ebreak:
		return e**(-index_lo)
	else:
		return e**(-index_hi) * emax**(index_hi - index_lo)

# Calculate photon flux from energy flux
def calc_flux(spectrum, e_flux):

	if spectrum['type'] == 'File':
		raise RuntimeError("If using spectrum .dat file, must define photon flux in input .yaml file.")
	elif spectrum['type'] == 'Mono':
		return e_flux * 6.24150647e8 / spectrum['energy']
	elif spectrum['type'] == 'Band':
		band_func = lambda e: band_function(e, spectrum['alpha'], spectrum['beta'], spectrum['ebreak'])
		return e_flux * 6.24150647e8 / integrate.quad(band_func, spectrum['energy_min'], spectrum['energy_max'])[0]
	elif spectrum['type'] == 'Comptonized':
		comp_func = lambda e: comp_function(e, spectrum['index'], spectrum['epeak'])
		return e_flux * 6.24150647e8 / integrate.quad(comp_func, spectrum['energy_min'], spectrum['energy_max'])[0]
	elif spectrum['type'] == 'PowerLaw':
		pl_func = lambda e: powerlaw_function(e, spectrum['index'])
		return e_flux * 6.24150647e8 / integrate.quad(pl_func, spectrum['energy_min'], spectrum['energy_max'])[0]
	elif spectrum['type'] == 'BrokenPowerLaw':
		bpl_func = lambda e: brokenpowerlaw_function(e, spectrum['ebreak'], spectrum['index_lo'], spectrum['index_hi'], spectrum['energy_max'])
		return e_flux * 6.24150647e8 / integrate.quad(bpl_func, spectrum['energy_min'], spectrum['energy_max'])[0]

# Defines zenith, azimuth, and flux of source
def define_angles_flux(inputs, spectrum):

	if inputs['coordinate_system'] == 'local':
		if 'zenith' in inputs:
			if type(inputs['zenith']) == int or type(inputs['zenith']) == float:
				zenith = float(inputs['zenith'])
			elif type(inputs['zenith']) == list:
				zenith = float(np.random.choice(inputs['zenith']))
			else:
				raise RuntimeError("'zenith' must be a number or list of numbers.")
		elif 'zenith_min' in inputs and 'zenith_max' in inputs:
			z_min = float(inputs['zenith_min'])
			z_max = float(inputs['zenith_max'])
			zenith = np.random.uniform(z_min, z_max)
		else:
			raise RuntimeError("Must specify zenith angle(s) in input .yaml file.")

		if 'azimuth' in inputs:
			if type(inputs['azimuth']) == int or type(inputs['azimuth']) == float:
				azimuth = float(inputs['azimuth'])
			elif type(inputs['azimuth']) == list:
				azimuth = float(np.random.choice(inputs['azimuth']))
			else:
				raise RuntimeError("'azimuth' must be a number or list of numbers.")
		elif 'azimuth_min' in inputs and 'azimuth_max' in inputs:
			a_min = float(inputs['azimuth_min'])
			a_max = float(inputs['azimuth_max'])
			azimuth = np.random.uniform(a_min, a_max)
		else:
			raise RuntimeError("Must specify azimuthal angle(s) in input .yaml file.")

	else:
		raise RuntimeError("Only detector coordinates are supported for now. 'coordinate_system' in input .yaml file must be 'local'.")

	if 'ph_flux' in inputs:
		if type(inputs['ph_flux']) == int or type(inputs['ph_flux']) == float:
			flux = float(inputs['ph_flux'])
		elif type(inputs['ph_flux']) == list:
			flux = float(np.random.choice(inputs['ph_flux']))
		else:
			raise RuntimeError("'ph_flux' must be a number or list of numbers.")
	elif 'ph_flux_min' in inputs and 'ph_flux_max' in inputs:
		flux_min = float(inputs['ph_flux_min'])
		flux_max = float(inputs['ph_flux_max'])
		flux = np.random.uniform(flux_min, flux_max)
	elif 'e_flux' in inputs:
		if type(inputs['e_flux']) == int or type(inputs['e_flux']) == float:
			flux = float(calc_flux(spectrum, inputs['e_flux']))
		elif type(inputs['e_flux']) == list:
			e_flux = float(np.random.choice(inputs['e_flux']))
			flux = float(calc_flux(spectrum, e_flux))
		else:
			raise RuntimeError("'e_flux' must be a number or list of numbers.")
	elif 'e_flux_min' in inputs and 'e_flux_max' in inputs:
		e_flux_min = float(inputs['e_flux_min'])
		e_flux_max = float(inputs['e_flux_max'])
		e_flux = np.random.uniform(e_flux_min, e_flux_max)
		flux = float(calc_flux(spectrum, e_flux))
	else:
		raise RuntimeError("Must specify flux(es) in input .yaml file.")

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

	return spectrum_text, lightcurve_text, spectrum

# Creates source file
def make_source_file(event, filename, geometry, shield_counts, z, a, flux, spectrum_text, lightcurve_text, coordsys):

	if coordsys == 'local':
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
	else:
		raise RuntimeError("Only detector coordinates are supported for now. 'coordinate_system' in input .yaml file must be 'local'.")


def main():
	input_file = parse_args().yaml
	inputs = read_yaml(input_file)
	input_path, output_path = define_paths(inputs)
	write_readme(inputs, input_path, output_path)
	events, lightcurves, spectra = make_event_dict(input_path, inputs['mix_or_match'], inputs['spectrum_type'])
	for event in events:
		source_filename = output_path + event + '.source'
		spectrum_text, lightcurve_text, spectrum = lightcurve_spectrum_text(event, lightcurves[event], spectra[event], input_path)
		zenith_angle, azimuthal_angle, flux = define_angles_flux(inputs, spectrum)
		make_source_file(event, source_filename, inputs['geometry_path'], inputs['shield_counts'], zenith_angle, azimuthal_angle, flux, spectrum_text, lightcurve_text, inputs['coordinate_system'])

if __name__ == "__main__":
    main()
