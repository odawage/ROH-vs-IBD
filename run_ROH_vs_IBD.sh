#!/bin/bash
#SBATCH --ntasks=2              # 1 core(CPU)
#SBATCH --nodes=1                # Use 1 node
#SBATCH --job-name=ROH_vs_IBD   # sensible name for the job
#SBATCH --mem=15G                 # Default memory per CPU is 3GB.
# SBATCH --partition=hugemem-avx2     # Use the smallmem-partition for jobs requiring < 10 GB RAM.
#SBATCH --array=0-39%20
## 0-39%
#SBATCH -o ./Reports/ROH_IBD_%A_%a.out    #Standar output message
#SBATCH -e ./Reports/ROH_IBD_%A_%a.err    #Standar error message

# If you would like to use more please adjust this.

# names=("NE200_100cm_rep0")

datasets=($(awk '{print $1}' Scenarios))
densities=($(awk '{print $2}' Scenarios))
parameters=($(awk '{print $3}' Scenarios))


dataset="${datasets[SLURM_ARRAY_TASK_ID]}"
dens="${densities[SLURM_ARRAY_TASK_ID]}"
parameter="${parameters[SLURM_ARRAY_TASK_ID]}"

echo "I am running" $dataset $dens $parameter



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

Rscript --vanilla ROH_vs_IBD.R "$dataset" "$dens" "$parameter"

