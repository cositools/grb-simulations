# Example Notebooks

This directory contains example Jupyter notebooks that download GBM data, generate source models from GBM data, produce input files for MEGAlib, simulate transients using MEGAlib, combine the simulations with background, and extract the ACS hits from the simulations.  

- `download_gbm_data`: Downloads data from GBM's Burst Catalog, as well as the TTE data
- `generate_source_models`: Produces lightcurves binned using Bayesian blocks from GBM TTE data, and creates spectral models from the GBM best spectral fits
- `create_megalib_inputs`: Generate MEGAlib input files using source models
- `run-simulations`: Simulate sources using MEGAlib
- `produce-acs-datasets`: Extract ACS hits from MEGAlib simulations, and combine with background