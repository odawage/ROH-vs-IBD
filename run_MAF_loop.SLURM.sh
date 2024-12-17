#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1                  # Number of tasks (cores)
#SBATCH --cpus-per-task=20
#SBATCH --job-name=ROH_MAF   # sensible name for the job
#SBATCH --mem=50G                 # Default memory per CPU is 3GB.
#SBATCH --partition=hugemem-avx2,orion,hugemem   # Partition name
#SBATCH --exclude=cn-14,cn-11
#SBATCH -o ROH_MAF.out    #Standard output message

# If you would like to use more please adjust this.



## Below you can put your scripts
# If you want to load module

module load Miniconda3
module load PLINK/1.9b_6.17-x86_64 
module load R
module list

## Prpare the environment for conda activation in the computing node:

#eval "$(conda shell.bash hook)"

##Load conda environment 

#conda activate PLINK2
#echo "I am working with the conda env: "$CONDA_PREFIX



## Run 

Rscript --vanilla MAF_loop.R
