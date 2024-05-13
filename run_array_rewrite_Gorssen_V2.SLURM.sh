#!/bin/bash
#SBATCH --ntasks=2              # 1 core(CPU)
#SBATCH --nodes=1                # Use 1 node
#SBATCH --job-name=Gorssen   # sensible name for the job
#SBATCH --mem=75G                 # Default memory per CPU is 3GB.
#SBATCH --partition=hugemem-avx2     # Use the smallmem-partition for jobs requiring < 10 GB RAM.
#SBATCH --exclude=cn-14
#SBATCH --array=0

#SBATCH -o ./Reports/Gorssen_%j_%a.out    #Standar output message
#SBATCH -e ./Reports/Gorssen_%j_%a.err    #Standar error message

# If you would like to use more please adjust this.

names=("NE200_100cm_rep0" "NE200_100cm_rep0_50" "NE200_100cm_rep1" "NE200_100cm_rep1_50" "NE200_100cm_rep2" "NE200_100cm_rep2_50" "NE100_100cm_rep1" "NE100_100cm_rep1_50" "NE100_100cm_rep2" "NE100_100cm_rep2_50")

# Get the name based on the array index
name="${names[SLURM_ARRAY_TASK_ID]}"



## Below you can put your scripts
# If you want to load module

module purge 
module load Miniconda3 
module load R
module load PLINK/1.9b_6.17-x86_64


module list

##Prpare the environment for conda activation in the computing node:

eval "$(conda shell.bash hook)"

##Load conda environment 

conda activate PLINK_R
echo "I am working with the conda env: "$CONDA_PREFIX



##Run 

Rscript --vanilla Array_rewrite_Gorssen_V2_vask.R "$name"
