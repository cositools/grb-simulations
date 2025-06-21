import os
import shutil
import glob
import logging

logger = logging.getLogger(__name__)

def run_cosima(file, output_dir, zipped=True, overwrite=False):
	'''
	Run cosima on a .source file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to .source file
	output_dir : pathlib.PosixPath
		Path to directory to store output .sim or .sim.gz file
	zipped : bool, optional
		Whether to zip .sim file
	overwrite : bool, optional
		Whether to overwrite existing .sim or .sim.gz file if it exists
	'''

	name = file.stem
	source_file = file.name

	if not overwrite:

		unzipped_sim_path = output_dir / name + '.sim'
		zipped_sim_path = output_dir / name + '.sim.gz'

		if unzipped_sim_path.exists() or zipped_sim_path.exists():

			logger.warning(f'Skipping {source_file} because .sim file already exists.')
			return

	output_dir.mkdir(parents=True, exist_ok=True)
	os.chdir(file.parent)

	if zipped:

		original_files = set(file.parent.glob(f'{name}*.sim.gz'))
		os.system(f'cosima {source_file} > {name}.log')

		new_files = set(file.parent.glob(f'{name}*.sim.gz')) - original_files
		for this_file in new_files:
			shutil.move(this_file, output_dir / name + '.sim.gz')

	else:

		original_files = set(file.parent.glob(f'{name}*.sim'))
		os.system(f'cosima -u {source_file} > {name}.log')

		new_files = set(file.parent.glob(f'{name}*.sim')) - original_files
		for this_file in new_files:
			shutil.move(this_file, output_dir / name + '.sim')

	if os.path.exists(file.parent / name + '.log'):
		os.remove(file.parent / name + '.log')

	if (file.parent / 'absorptions').exists() and (file.parent / 'absorptions').is_dir():
		shutil.rmtree(file.parent / 'absorptions')

def revan(file, output_dir, mass_model, config_file, overwrite=False):
	'''
	Run revan on a .sim or .sim.gz file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to .sim or .sim.gz file
	output_dir : pathlib.PosixPath
		Path to directory to store output .tra.gz file
	mass_model : pathlib.PosixPath
		Path to analysis mass model
	config_file : pathlib.PosixPath
		Path to revan configuration file
	overwrite : bool, optional
		Whether to overwrite existing .tra.gz file if it exists
	'''

	name = file.stem
	sim_file = file.name

	if not overwrite:

		tra_path = output_dir / name + '.tra.gz'

		if tra_path.exists():

			logger.warning(f'Skipping {sim_file} because .tra.gz file already exists.')
			return

	output_dir.mkdir(parents=True, exist_ok=True)
	os.chdir(file.parent)

	original_files = set(file.parent.glob(f'{name}*.tra.gz'))
	os.system(f'revan -f {sim_file} -g {str(mass_model)} -c {str(config_file)} -n -a > {name}.log')

	new_files = set(file.parent.glob(f'{name}*.tra.gz')) - original_files
	for this_file in new_files:
		shutil.move(this_file, output_dir / name + '.tra.gz')

	if os.path.exists(file.parent / name + '.log'):
		os.remove(file.parent / name + '.log')

	if (file.parent / 'absorptions').exists() and (file.parent / 'absorptions').is_dir():
		shutil.rmtree(file.parent / 'absorptions')

def mimrec(file, output_dir, mass_model, config_file, overwrite=False):
	'''
	Run mimrec on a .tra.gz file.

	Parameters
	----------
	file : pathlib.PosixPath
		Path to .tra.gz file
	output_dir : pathlib.PosixPath
		Path to directory to store output .tra file
	mass_model : pathlib.PosixPath
		Path to analysis mass model
	config_file : pathlib.PosixPath
		Path to mimrec configuration file
	overwrite : bool, optional
		Whether to overwrite existing .tra file if it exists
	'''

	name = file.stem
	tra_file = file.name

	if not overwrite:

		extracted_path = output_dir / name + '.extracted.tra'

		if extracted_path.exists():

			logger.warning(f'Skipping {tra_file} because .extracted.tra file already exists.')
			return

	output_dir.mkdir(parents=True, exist_ok=True)
	os.chdir(file.parent)

	original_files = set(file.parent.glob(f'{name}*.tra'))
	os.system(f'mimrec -f {tra_file} -g {str(mass_model)} -c {str(config_file)} -n -x > {name}.log')

	new_files = set(file.parent.glob(f'{name}*.tra')) - original_files
	for this_file in new_files:
		shutil.move(this_file, output_dir / name + '.extracted.tra')

	if os.path.exists(file.parent / name + '.log'):
		os.remove(file.parent / name + '.log')

	if (file.parent / 'absorptions').exists() and (file.parent / 'absorptions').is_dir():
		shutil.rmtree(file.parent / 'absorptions')

def simulate(source_dir, output_dir, zipped=True, overwrite=False):
	'''
	Run cosima on all .source files in a directory.

	Parameters
	----------
	source_dir : pathlib.PosixPath
		Path to directory housing .source files
	output_dir : pathlib.PosixPath
		Path to directory to store output .sim or .sim.gz files
	zipped : bool, optional
		Whether to zip .sim files
	overwrite : bool, optional
		Whether to overwrite existing .sim or .sim.gz files if they exist
	'''

	for source_file in source_dir.iterdir():

		if source_file.suffix == '.source':

			run_cosima(source_file, output_dir, zipped, overwrite)

def reconstruct(sim_dir, output_dir, mass_model, config_file, overwrite=False):
	'''
	Run revan on all .sim and .sim.gz files in a directory.

	Parameters
	----------
	sim_dir : pathlib.PosixPath
		Path to directory housing .sim and/or .sim.gz files
	output_dir : pathlib.PosixPath
		Path to directory to store output .tra.gz files
	mass_model : pathlib.PosixPath
		Path to analysis mass model
	config_file : pathlib.PosixPath
		Path to revan configuration file
	overwrite : bool, optional
		Whether to overwrite existing .tra.gz files if they exist
	'''

	for sim_file in sim_dir.iterdir():

		if sim_file.suffix == '.sim' or sim_file.suffix == '.sim.gz':

			run_revan(sim_file, output_dir, mass_model, config_file, overwrite)

def extract(tra_dir, output_dir, mass_model, config_file, overwrite=False):
	'''
	Run mimrec on all .tra.gz files in a directory.

	Parameters
	----------
	tra_dir : pathlib.PosixPath
		Path to directory housing .tra.gz files
	output_dir : pathlib.PosixPath
		Path to directory to store output extracted .tra files
	mass_model : pathlib.PosixPath
		Path to analysis mass model
	config_file : pathlib.PosixPath
		Path to mimrec configuration file
	overwrite : bool, optional
		Whether to overwrite existing .tra files if they exist
	'''

	for tra_file in tra_dir.iterdir():

		if tra_file.suffix == '.tra.gz':

			run_mimrec(tra_file, output_dir, mass_model, config_file, overwrite)
