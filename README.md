# COSI GRB Group Simulations

This repository contains code to run COSI simulations of GRBs and other short-duration transients in bulk. 

Simulations for the GRB group run using this pipeline should be uploaded [here](https://drive.google.com/drive/u/0/folders/11_qUIzQx3oGTjrb6voim0GB_EgXny9co).  

If you have questions or encounter any problems, please open an issue.  

## Installation Instructions     
1. Clone the repository. From the command line, run:  
`git clone https://github.com/cositools/grb-simulations.git`  
2. Navigate into the project directory.  
`cd grb-simulations`  
3. Create and activate a conda environment (this is optional, but highly recommended).  
<pre> ```conda create -n cosi-grb-sims python=3.11
conda activate cosi-grb-sims``` </pre>  
4. Install the cosiburstpy package, ideally within a conda environment.  
`pip install .`  