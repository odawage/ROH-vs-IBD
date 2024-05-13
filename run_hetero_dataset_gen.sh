#!/bin/bash
#SBATCH --ntasks=20              # 1 core(CPU)
#SBATCH --nodes=1                # Use 1 node
#SBATCH --job-name=Froh_Het   # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB.
# SBATCH --partition=hugemem-avx2     # Use the smallmem-partition for jobs requiring < 10 GB RAM
#SBATCH -o ./Reports/Froh_Het_%j.out    #Standar output message
#SBATCH -e ./Reports/Froh_Het_%j.err    #Standar error message

# If you would like to use more please adjust this.

# names=("NE200_100cm_rep0")


# load modules
module load Miniconda3 
module list

## Prpare the environment for conda activation in the computing node:

eval "$(conda shell.bash hook)"

##Load conda environment 

conda activate PLINK2
echo "I am working with the conda env: "$CONDA_PREFIX



## Run 

Rscript --vanilla Hetero_dataset_gen_vask.R
