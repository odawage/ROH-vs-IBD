#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1                  # Number of tasks (cores)
#SBATCH --cpus-per-task=40
#SBATCH --job-name=ROH_err   # sensible name for the job
#SBATCH --mem=50G                 # Default memory per CPU is 3GB.
#SBATCH --partition=hugemem-avx2,orion,hugemem   # Partition name
#SBATCH --exclude=cn-14,cn-11
#SBATCH -o ROH_err.out    #Standard output message

# If you would like to use more please adjust this.

#names=("NE200_100cm_rep0" "NE200_100cm_rep1" "NE200_100cm_rep2" "NE100_100cm_rep1" "NE100_100cm_rep2")

# Get the name based on the array index
#name="${names[SLURM_ARRAY_TASK_ID]}"

#echo $name

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

Rscript --vanilla gen_err_loop.R
