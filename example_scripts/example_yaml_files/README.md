# Example .yaml and .dat files

This directory contains examples of the .yaml and .dat files needed to run the python scripts.

`example_download.yaml` is a sample .yaml file used as input for `download_gbm_data.py`.    
Parameters to be defined in the .yaml file:    
- output_path: path to directory to store downloaded files (e.g. 'gbm_data/')
- filters: parameters used to filter bursts (must be columns in [table](https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dfermigbrst&Action=More+Options))
- download: columns of GBM burst data to download (must be columns in [table](https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dfermigbrst&Action=More+Options), must include trigger_name and bcat_detector_mask)    

`example_gbm_to_megalib.yaml` is a sample .yaml file used as input for `gbm_to_megalib.py`.      
Parameters to be defined in the .yaml file:      
- input_path: path to directory housing downloaded GBM data (e.g. 'gbm_data/')      
- output_path: path to directory to store MEGAlib source inputs (e.g. 'MEGAlib_source_inputs/')       
- plot_path: path to directory to store lightcurve plots      
- source_time_range: time range of source lightcurve (s)     
- background_time_range: time range to use to fit background (s)     
- nai_energy_range: energy range of GBM's NaI detectors used for lightcurve     
- bgo_energy_range: energy range of GBM's BGO detectors used for lightcurve     
- cosi_energy_range: energy range for COSI MEGAlib simulations 
- min_duration: minimum duration of events to reject once binned using Bayesian blocks (s), should match maximum expected event duration      

`example_source_file.yaml` is a sample .yaml file used as input for `sample_source_file_generation.py`. When `sample_source_file_generation.py` is run, it takes spectral (.yaml & .dat) and lightcurve (.dat) files from the directory specified in the input .yaml file (`example.yaml`) and writes the spectral info from the spectrum .yaml file (if spectrum_type: 'yaml') or the path to the .dat file into the MEGAlib .source file (spectrum_type: 'dat'), as well as the path to the lightcurve .dat file.      
Parameters to be defined in the .yaml file:  
- input_path: path to directory housing source inputs (e.g. 'MEGAlib_source_inputs/')   
- output_path: path to directory to store .source files (e.g. 'MEGAlib_source_files/')  
- source_input_path: path to directory where .source files should look for lightcurve and spectral .dat files (optional, only needed if different than input path)      
- mass_model_path: path to mass model  
- spectrum_type: file type to use for spectrum ('yaml' to use built in MEGAlib spectral models with parameters defined in .yaml files or 'dat' to use .dat files with spectral shape)  
- mix_or_match: whether to mix or match lightcurves and spectra ('match' will use lightcurves and spectra from same event, 'mix' will randomly match lightcurves with spectra from input directory)  
- zenith: incidence (zenith) angle or list of incidence (zenith) angles (between 0 and 180 where 90 degrees is on-axis)    
- zenith_min: minimum incidence (zenith) angle (between 0 and 180 where 90 degrees is on-axis, only used if 'zenith' is not defined)     
- zenith_max: maximum incidence (zenith) angle (between 0 and 180 where 90 degrees is on-axis, only used if 'zenith' is not defined)     
- azimuth: azimuthal angle or list of azimuthal angles (between 0 and 360)    
- azimuth_min: minimum azimuthal angle (between 0 and 360, only used if 'azimuth' is not defined)     
- azimuth_max: maximum azimuthal angle (between 0 and 360, only used if 'azimuth' is not defined)     
- ph_flux: flux or list of fluxes in ph/cm<sup>2</sup>/s   
- ph_flux_min: minimum flux in ph/cm<sup>2</sup>/s (only used if 'ph_flux' is not defined)    
- ph_flux_max: maximum flux in ph/cm<sup>2</sup>/s (only used if 'ph_flux' is not defined)    
- e_flux: flux or list of fluxes in erg/cm<sup>2</sup>/s (only used if 'ph_flux', 'ph_flux_min', & 'ph_flux_max' are not defined)      
- e_flux_min: minimum flux in erg/cm<sup>2</sup>/s (only used if 'e_flux', 'ph_flux', 'ph_flux_min', & 'ph_flux_max' are not defined)    
- e_flux_max: maximum flux in erg/cm<sup>2</sup>/s (only used if 'e_flux', 'ph_flux', 'ph_flux_min', & 'ph_flux_max' are not defined)    
- shield_counts: whether to include line in .source files to store background shield counts in MEGAlib .sim files (optional, 'y' or True to include)   
- coordinate_system: Coordinate system of .source files ('local' for detector coordinates or 'galactic' for galactic coordinates), galactic is not yet supported   

The `example-spectraltype_spectrum.yaml` files are examples of each of the MEGAlib spectral models supported by `sample_source_file_generation.py`: monoenergetic, Band function, Comptonized, power law, and broken power law. The parameters needed depend on the model. These files must end in '_spectrum.yaml' and be located in the `input_path/` directory indicated in the input .yaml file when running `source_file_generation.py`.      

`example_spectrum.dat` is a sample spectral .dat file which can be used in place of the .yaml files. This allows for custom spectral shapes to be used instead of the built in MEGAlib functions. Only the shape is important, not the amplitude, since the flux is defined in the .source file. These files must end in '_spectrum.dat' and be located in the `input_path/` directory indicated in the input .yaml file when running `source_file_generation.py`.      

`example_lightcurve.dat` is a sample lightcurve file. Only the shape is important, not the amplitude, since the flux is defined in the .source file. These files must end in '_lightcurve.dat' and be located in the `input_path/` directory indicated in the input .yaml file when running `source_file_generation.py`.    

`example_cosima.yaml` is a sample .yaml file used as input for `run_cosima.py`.     
Parameters to be defined in the .yaml file:   
- input_path: path to directory housing MEGAlib .source files (e.g. 'MEGAlib_source_files/')      
- output_path: path to directory to store .sim files (e.g. 'MEGAlib_outputs/')     

`example_revan.yaml` is a sample .yaml file used as input for `run_revan.py`.     
Parameters to be defined in the .yaml file:   
- input_path: path to directory housing .sim files (e.g. 'MEGAlib_outputs/')      
- output_path: path to directory to store .tra files (e.g. 'MEGAlib_outputs/')    
- mass_model_path: path to instrument mass model (e.g. 'COSISMEX.O64.geo.setup')      
- config_file_path: path to revan configuration file (e.g. 'SMEXv12.Continuum.HEALPixO3.binnedimaging.revan.cfg')      

`example_mimrec.yaml` is a sample .yaml file used as input for `run_mimrec.py`.     
Parameters to be defined in the .yaml file:   
- input_path: path to directory housing .tra files (e.g. 'MEGAlib_outputs/')      
- output_path: path to directory to store extracted .tra files (e.g. 'MEGAlib_outputs/')    
- mass_model_path: path to instrument mass model (e.g. 'COSISMEX.O64.geo.setup')      
- config_file_path: path to mimrec configuration file (e.g. 'SMEXv12.Continuum.HEALPixO3.binnedimaging.mimrec.cfg')     

`example_trigger_algorithm.yaml` is a sample .yaml file used as input for `trigger_algorithm_list_generation.py`.    
Parameters to be defined in the .yaml file:    
- source_path: path to directory housing source .sim files (e.g. 'MEGAlib_outputs/')  
- output_path: path to directory to store output (e.g. 'trigger_inputs/')    
- source_file_path: path to directory housing source files used to simulate .sim files (e.g. 'MEGAlib_source_files/')
- mass_model_path: path to geometry file, must match the geometry used to create simulations    
- background_type: whether to randomly select background regions for each source ('random') or use same background file for all sources ('file')  
- background_path: if background_type is 'random', path to directory housing background .sim files, and if background_type is 'file', the path to the concatenated background file (e.g. 'MEGAlib_backgrounds/')    
- background_number: number of background files for each component (only necessary if background_type is 'random')    
- background_file_length: length in s of each background file (only necessary if background_type is 'random')    
- background_time: amount of background in s to add to each source, source will be added to the middle of the background time interval    
- background_components: list of background components to include, file names must begin with component names (only necessary if background_type is 'random')    
- background_file_type: whether the background files begin one after another in time ('sequential') or if they begin at the same time ('simultaneous') (only necessary if background_type is 'random')    
- mass_model_version: version of mass model   
