# COSI GRB Group Simulations

This repository contains the code and inputs for COSI GRB simulations, as well as other short-duration transients which will need to be classified and investigated by the GRB group.

## Steps to Running Simulations      
1. Create source input files, either using data observed by other telescopes (e.g. Fermi-GBM) or theoretical models. `download_gbm_data.py` is an example script that can be used to download lightcurves and spectra from Fermi-GBM, and `gbm_to_megalib.py` may be used to convert the downloaded lightcurves and spectra to the formats required to generate the MEGAlib .source files. These files should be uploaded to `MEGAlib_source_inputs/`.       
2. Generate MEGAlib .source files for transient events. `source_file_generation.py` may be used to generate these files. They are too large to store on GitHub and should be uploaded [here](https://drive.google.com/drive/folders/1SAxRKPWObcZQ8T93njdKWMTDuI4_-Bx3).      
3. Simulate the .source files using cosima, reconstruct the resulting .sim files using revan, and perform event selections on the resulting .tra files using mimrec. `run_cosima.py`, `run_revan.py`, and `run_mimrec.py` are example scripts that can be used to run cosima, revan, and mimrec, respectively. For onboard trigger algorithm simulations or anti-coincidence shield (BGO) only simulations, only the cosima step is necessary. The MEGAlib simulations are too large to store on GitHub and should be uploaded [here](https://drive.google.com/drive/folders/1aPKtB7uTr5LNG3xGbeCpIPTbAW6ibvk3).   
4. Combine source simulations with background simulations. Depending on the analysis, this may be done at either the .sim or .tra file level. There are MEGAlib tools (e.g. mtraconcatter) which do this (see the [MEGAlib documentation](https://megalibtoolkit.com/documentation.html)), and if creating event lists to be used to test the onboard trigger algorithm, this is done in `trigger_algorithm_list_generation.py`.

## Scripts     
- `download_gbm_data.py` downloads Fermi-GBM spectra and lightcurves given an input .yaml file (an example can be found in `examples/`) and stores the data in a specified output directory.   
To run from the command line: `python download_gbm_data.py -y input-file.yaml`.    
- `gbm_to_megalib.py` creates MEGAlib inputs (spectral .yaml files and lightcurve .dat files binned using Bayesian blocks) from GBM data given an input .yaml file (an example can be found in `examples/`). Currently, this only produces MEGAlib input files for events with Band function or Comptonized best fit spectra.      
To run from the command line: `python gbm_to_megalib.py -y input-file.yaml`.    
- `source_file_generation.py` creates source file samples based an on an input .yaml file (an example can be found in `examples/`). This draws lightcurves and spectra from the samples in `MEGAlib_source_inputs/` and create .source files in a specified directory. The .source files are created with the defined mass model path, fluxes, and incidence angles. Currently, this only produces .source files in detector coordinates, but the ability to create .source files in Galactic coordinates with orientation files will be added later, as well as polarization.          
To run from the command line: `python source_file_generation.py -y input-file.yaml`.      
- `run_cosima.py` runs cosima on all of the .source files in a specified directory and produces simulation (.sim) files. You must be on the main branch of MEGAlib to run cosima.       
To run from the command line: `python run_cosima.py -y input-file.yaml`.  
- `run_revan.py` runs revan on all of the .sim files in a specified directory and produces reconstructed simulation (.tra) files. You must be on the dee2022 branch of MEGAlib to run revan.    
To run from the command line: `python run_revan.py -y input-file.yaml`.    
- `run_mimrec.py` runs mimrec on all of the .tra files in a specified directory and produces reconstructed simulation (.tra) files with event selections. You must be on the dee2022 branch of MEGAlib to run mimrec.     
To run from the command line: `python run_mimrec.py -y input-file.yaml`.    
- `trigger_algorithm_list_generation.py` creates event lists with the time and energy of each hit in the BGO and GeDs based on an input .yaml file (an example can be found in `examples/`). This draws source and background simulations (.sim files) from specified input directories and creates trigger algroithm files in a specified output directory.         
To run from the command line: `python trigger_algorithm_list_generation.py -y input-file.yaml`.      

## Directories          
- `examples/` contains sample .yaml files needed for `download_gbm_data.py`, `gbm_to_megalib.py`, `source_file_generation.py`, and `trigger_algorithm_list_generation.py`, as well as sample spectral and lightcurve files. See the README file in `examples/` for more details.       
- `MEGAlib_source_inputs/` contains the information necessary to generate a given simulation set in COSI. This is sub-divided into folders by event type, currently including gamma-ray bursts (grbs), magnetar giant flares (mgfs), soft gamma-ray repeater short bursts (sgr bursts), solar flares (solar_flares), and terrestrial gamma-ray flashses (tgfs). Each of these is further subdivided into the specific input samples generated by various team members. These will generally be spectral (.dat/.yaml) and lightcurve (.dat) files, as well as relevant images (.pdf, .png, etc.). Each of these individual samples should contain a text or catalog file describing the sample as well as the responsible team member. For use with `source_file_generation.py`, spectral files must be in the form `event-name_spectrum.dat` or `event-name_spectrum.yaml` and lightcurves must be in the form `event-name_lightcurve.dat`. Examples of spectral and lightcurve .yaml files can be found in `examples/`, and information about the spectral .dat file format can be found in the [MEGAlib documentation](https://megalibtoolkit.com/documentation.html).       

## Uploaded Files
- MEGAlib source files (outputs of `sample_source_file_generation.py`) can be found [here](https://drive.google.com/drive/folders/1SAxRKPWObcZQ8T93njdKWMTDuI4_-Bx3).     
- Source simulations of GRBs and other transients (outputs of `run_cosima.py`, `run_revan.py`, and `run_mimrec.py`) can be found [here](https://drive.google.com/drive/folders/1aPKtB7uTr5LNG3xGbeCpIPTbAW6ibvk3), and background simulations can be found [here](https://drive.google.com/drive/folders/1OTN-_8gUxedueEbL3mPeh0_0kh7e9kKF).      
- Event list files used to test the onboard trigger algorithm (outputs of `trigger_algorithm_list_generation.py`) can be found [here](https://drive.google.com/drive/folders/1ibXOrTFD2-ntYy54HbWRXj7YYmueEl0e).     
