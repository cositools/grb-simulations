# COSI GRB Group Simulations

This repository contains the code and inputs for COSI GRB simulations, as well as other short-duration transients which will need to be classified and investigated by the GRB group. This is a work in progress; please open an issue if you notice any problems/mistakes or missing features, and contact Eliza (eneights@gwu.edu) with any questions.

Simulations run by Eliza using this pipeline can be found [here](https://drive.google.com/drive/u/0/folders/11_qUIzQx3oGTjrb6voim0GB_EgXny9co).

## Steps to Running Simulations      
1. Create source input files, either using data observed by other telescopes (e.g. Fermi-GBM) or theoretical models. `download_gbm_data.py` is an example script that can be used to download lightcurves and spectra from Fermi-GBM.  `gbm_to_megalib.py` may be used to convert the downloaded lightcurves and spectra to the formats required to generate the MEGAlib .source files.         
2. Generate MEGAlib .source files for transient events. `source_file_generation.py` may be used to generate these files. For use with `source_file_generation.py`, spectral files must be in the form `event-name_spectrum.dat` or `event-name_spectrum.yaml` and lightcurves must be in the form `event-name_lightcurve.dat`. Examples of spectral and lightcurve .yaml files can be found in `example_yaml_files/`, and information about the spectral .dat file format can be found in the [MEGAlib documentation](https://megalibtoolkit.com/documentation.html).      
3. Simulate the .source files using cosima, reconstruct the resulting .sim files using revan, and perform event selections on the resulting .tra files using mimrec. `run_cosima.py`, `run_revan.py`, and `run_mimrec.py` are example scripts that can be used to run cosima, revan, and mimrec, respectively. For onboard trigger algorithm simulations or anti-coincidence shield (BGO) only simulations, only the cosima step is necessary.    
4. Combine source simulations with background simulations. Depending on the analysis, this may be done at either the .sim or .tra file level. There are MEGAlib tools (e.g. mtraconcatter) which do this (see the [MEGAlib documentation](https://megalibtoolkit.com/documentation.html)), and `trigger_algorithm_list_generation.py` can be used to create input files for an implementation of the onboard trigger algorithm. 

## Directories          
- `example_scripts/` contains scripts that can be used to download Fermi-GBM spectra and lightcurves, bin the lightcurves using Bayesian blocks and convert the spectra and lightcurves into MEGAlib input files, create MEGAlib .source files, run simulations in MEGAlib, and create input files for the onboard trigger algorithm.              
- `grb_simulator/` contains the code needed to run simulations of transients in bulk.           

