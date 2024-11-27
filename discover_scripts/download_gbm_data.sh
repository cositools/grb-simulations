#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH -o output.%j
#SBATCH -e error.%j
#SBATCH --account=j1042
#SBATCH --job-name=download-sgrbs
#SBATCH --ntasks=1
#SBATCH --workdir=$HOME/slurm_scripts/download_gbm_data/

module add anaconda/py3.9
conda activate gbm

python download_gbm_data.py -y download_gbm_sgrbs.yaml