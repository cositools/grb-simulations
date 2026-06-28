# COSI GRB Group Simulations

This repository contains code to run population-level COSI simulations of GRBs and other short-duration transients. 

Simulations for the GRB group run using this pipeline should be uploaded [here](https://drive.google.com/drive/u/0/folders/11_qUIzQx3oGTjrb6voim0GB_EgXny9co).  

If you have questions or encounter any problems, please open an issue.  

## Installation Instructions     
1. Clone the repository. From the command line, run:  
`git clone https://github.com/cositools/grb-simulations.git`  
2. Navigate into the project directory.  
`cd grb-simulations`  
3. Optional but highly recommended: Create and activate a conda environment.  
`conda create -n cosi-grb-sims python=3.12`  
`conda activate cosi-grb-sims`  
4. Install the cosiburstpy package.  
   a. If you do not need to work with Fermi-GBM data, run:  
      `pip install .`  
   b. If you are using Fermi-GBM data, run:  
       `pip install -e ".[fermi-gbm]"`  
       `cosiburstpy setup"`
