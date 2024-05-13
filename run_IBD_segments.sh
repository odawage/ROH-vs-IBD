#!/bin/bash
#SBATCH --ntasks=1              # 1 core(CPU)
#SBATCH --nodes=1                # Use 1 node
#SBATCH --job-name=IBD_seg   # sensible name for the job
#SBATCH --mem=15G                 # Default memory per CPU is 3GB.
#SBATCH --partition=orion              #hugemem-avx2      Use the smallmem-partition for jobs requiring < 10 GB RAM.
#SBATCH --array=0-4

#SBATCH -o ./Reports/IBD_seg_%A_%a.out    #Standar output message
#SBATCH -e ./Reports/IBD_seg_%A_%a.err    #Standar error message

# If you would like to use more please adjust this.

names=("NE200_100cm_rep0" "NE200_100cm_rep1" "NE200_100cm_rep2" "NE100_100cm_rep1" "NE100_100cm_rep2")

# Get the name based on the array index
name="${names[SLURM_ARRAY_TASK_ID]}"

echo $name

## Below you can put your scripts
# If you want to load module

module load Miniconda3 
module list

## Prpare the environment for conda activation in the computing node:

eval "$(conda shell.bash hook)"

##Load conda environment 

conda activate PLINK2
echo "I am working with the conda env: "$CONDA_PREFIX



## Run 

Rscript --vanilla IBD_segments.R "$name"
