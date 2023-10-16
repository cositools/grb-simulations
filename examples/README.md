# Example .yaml and .dat files

This directory contains examples of the .yaml and .dat files needed to run `sample_source_file_generation.py`.

`example.yaml` is a sample .yaml file used as input for `sample_source_file_generation.py`. 

The `example-spectraltype_spectrum.yaml` files are examples of each of the MEGAlib spectral models supported by `sample_source_file_generation.py`: monoenergetic, Band function, Comptonized, power law, and broken power law. These will need to be located in the `input_path/event_type/event_subtype` (or `input_path/event_type` if no event subtypes) directory indicated in the input .yaml file.

`example_spectrum.dat` is a sample spectral .dat file which can be used in place of the .yaml files. This allows for custom spectral shapes to be used instead of the built in MEGAlib functions. Only the shape is important, not the amplitude, since the flux is defined in the .source file.

`example_lightcurve.dat` is a sample lightcurve file. Only the shape is important, not the amplitude, since the flux is defined in the .source file.


