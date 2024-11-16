# Example .yaml and .dat files

This directory contains examples of the .yaml and .dat files needed to run the python scripts.

`example_download.yaml` is a sample .yaml file used as input for `download_gbm_data.py`.    
Parameters to be defined in the .yaml file:    
- paths:      
    - output: path to directory to store downloaded files (e.g. 'gbm_data/')       
- filters: parameters used to filter bursts (must be columns in [GBM burst table](https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dfermigbrst&Action=More+Options))      
- download: columns of GBM burst data to download (must be columns in [GBM burst table](https://heasarc.gsfc.nasa.gov/db-perl/W3Browse/w3table.pl?tablehead=name%3Dfermigbrst&Action=More+Options), must include trigger_name, bcat_detector_mask, and t90)     

`example_gbm_to_megalib.yaml` is a sample .yaml file used as input for `gbm_to_megalib.py`.      
Parameters to be defined in the .yaml file:      
- paths:    
    - input: path to directory housing downloaded GBM data (e.g. 'gbm_data/')      
    - output: path to directory to store MEGAlib source inputs (e.g. 'MEGAlib_source_inputs/')       
    - plots: path to directory to store lightcurve plots      
- time:     
    - source_range: time range of source lightcurve (s) with T90 beginning at t=0 s     
    - background_range: time range to use to fit background (s) with T90 beginning at t=0 s         
- energy:     
    - nai_range: energy range of GBM's NaI detectors used for lightcurve     
    - bgo_range: energy range of GBM's BGO detectors used for lightcurve     
    - cosi_range: energy range for COSI MEGAlib simulations      

`example_source_file.yaml` is a sample .yaml file used as input for `sample_source_file_generation.py`. When `sample_source_file_generation.py` is run, it takes spectral (.yaml & .dat) and lightcurve (.dat) files from the directory specified in the input .yaml file (`example.yaml`) and writes the spectral info from the spectrum .yaml file (if spectrum_type: 'yaml') or the path to the .dat file into the MEGAlib .source file (spectrum_type: 'dat'), as well as the path to the lightcurve .dat file.      
Parameters to be defined in the .yaml file:  
- paths:     
    - input: path to directory housing source inputs (e.g. 'MEGAlib_source_inputs/')   
    - output: path to directory to store .source files (e.g. 'MEGAlib_source_files/')  
    - source_input: path to directory where .source files should look for lightcurve and spectral .dat files (optional, only needed if different than input path)      
    - mass_model: path to mass model  
    - orientation: path to orientation file (optional, only needed if coordinate system is galactic)
- position:        
    - zenith: incidence (zenith) angle or list of incidence (zenith) angles (between 0 and 180 where 90 degrees is on-axis)    
    - zenith_range: list of minimum and maximum incidence (zenith) angles (between 0 and 180 where 90 degrees is on-axis, only used if 'zenith' is not defined)     
    - azimuth: azimuthal angle or list of azimuthal angles (between 0 and 360)    
    - azimuth_range: minimum and maximum azimuthal angles (between 0 and 360, only used if 'azimuth' is not defined)     
    - latitude: galactic latitude or list of galactic latitudes (between -90 and 90)       
    - longitude: galactic longitude or list of galactic longitudes (between 0 and 360)             
- spectra_and_lightcurves:       
    - photon_flux: flux or list of fluxes in ph/cm<sup>2</sup>/s   
    - photon_flux_range: list of minimum and maximum fluxes in ph/cm<sup>2</sup>/s (only used if 'photon_flux' is not defined)    
    - energy_flux: flux or list of fluxes in erg/cm<sup>2</sup>/s (only used if 'photon_flux' & 'photon_flux_range' are not defined)      
    - energy_flux_range: list of minimum and maximum fluxes in erg/cm<sup>2</sup>/s (only used if 'energy_flux', 'photon_flux', & 'photon_flux_range' are not defined)     
    - spectrum_type: file type to use for spectrum ('yaml' to use built in MEGAlib spectral models with parameters defined in .yaml files or 'dat' to use .dat files with spectral shape)  
    - mix_or_match: whether to mix or match lightcurves and spectra ('match' will use lightcurves and spectra from same event, 'mix' will randomly match lightcurves with spectra from input directory)       
    - start_time: time or list with minimum and maximum times in s at which source lightcurves will begin (optional, if not provided, will use the unmodified lightcurves for local coordinates, and choose random times within the orientation file for galactic coordinates)         
- general:       
    - shield_counts: whether to include line in .source files to store background shield counts in MEGAlib .sim files (optional, 'y' or True to include)   
    - coordinate_system: Coordinate system of .source files ('local' for detector coordinates or 'galactic' for galactic coordinates)   

The `example-spectraltype_spectrum.yaml` files are examples of each of the MEGAlib spectral models supported by `sample_source_file_generation.py`: monoenergetic, Band function, Comptonized, power law, and broken power law. The parameters needed depend on the model. These files must end in '_spectrum.yaml' and be located in the `paths: input` directory indicated in the input .yaml file when running `source_file_generation.py`.      

`example_spectrum.dat` is a sample spectral .dat file which can be used in place of the .yaml files. This allows for custom spectral shapes to be used instead of the built in MEGAlib functions. Only the shape is important, not the amplitude, since the flux is defined in the .source file. These files must end in '_spectrum.dat' and be located in the `paths: input` directory indicated in the input .yaml file when running `source_file_generation.py`.      

`example_lightcurve.dat` is a sample lightcurve file. Only the shape is important, not the amplitude, since the flux is defined in the .source file. These files must end in '_lightcurve.dat' and be located in the `paths: input` directory indicated in the input .yaml file when running `source_file_generation.py`.    

`example_cosima.yaml` is a sample .yaml file used as input for `run_cosima.py`.     
Parameters to be defined in the .yaml file:   
- paths:         
    - source_files: path to directory housing MEGAlib .source files (e.g. 'MEGAlib_source_files/')      
    - sim_files: path to directory to store .sim files (e.g. 'MEGAlib_outputs/')     
- cosima:          
    - zip: whether to gzip .sim files (optional, False to not zip files)
    - parallel: whether to run simulations in parallel (optional, True to run in parallel)
    - instances: number of mcosima instances to run (optional, only needed if 'parallel' is True)

`example_revan.yaml` is a sample .yaml file used as input for `run_revan.py`.     
Parameters to be defined in the .yaml file:   
- paths:         
    - sim_files: path to directory housing .sim files (e.g. 'MEGAlib_outputs/')      
    - tra_files: path to directory to store .tra files (e.g. 'MEGAlib_outputs/')    
    - mass_model: path to instrument mass model (e.g. 'COSISMEX.O64.geo.setup')      
- revan:          
    - config: path to revan configuration file (e.g. 'SMEXv12.Continuum.HEALPixO3.binnedimaging.revan.cfg')      

`example_mimrec.yaml` is a sample .yaml file used as input for `run_mimrec.py`.     
Parameters to be defined in the .yaml file:   
- paths:         
    - tra_files: path to directory housing .tra files (e.g. 'MEGAlib_outputs/')      
    - extracted_tra_files: path to directory to store extracted .tra files (e.g. 'MEGAlib_outputs/')    
    - mass_model: path to instrument mass model (e.g. 'COSISMEX.O64.geo.setup')      
- mimrec:          
    - config: path to mimrec configuration file (e.g. 'SMEXv12.Continuum.HEALPixO3.binnedimaging.mimrec.cfg')     

`example_trigger_algorithm.yaml` is a sample .yaml file used as input for `trigger_algorithm_list_generation.py`.    
Parameters to be defined in the .yaml file:    
- paths:       
    - input: path to directory housing source .sim files (e.g. 'MEGAlib_outputs/')  
    - output: path to directory to store output (e.g. 'trigger_inputs/')    
    - source_files: path to directory housing source files used to simulate .sim files (e.g. 'MEGAlib_source_files/')
    - mass_model: path to geometry file, must match the geometry used to create simulations but must be altered to not veto BGO hits
    - background: if background_type is 'random', path to directory housing background .sim files, and if background_type is 'file', the path to the concatenated background file (e.g. 'MEGAlib_backgrounds/')        
- background:       
    - type: whether to randomly select background regions for each source ('random') or use same background file for all sources ('file')  
    - number: number of background files for each component (only necessary if background_type is 'random')    
    - file_length: length in s of each background file (only necessary if background_type is 'random')    
    - time: amount of background in s to add to each source, source will be added to the middle of the background time interval    
    - components: list of background components to include, file names must begin with component names (only necessary if background_type is 'random')    
    - file_type: whether the background files begin one after another in time ('sequential') or if they begin at the same time ('simultaneous') (only necessary if background_type is 'random')    
- mass_model_version: version of mass model   
