#!/bin/bash
#SBATCH --ntasks=1              # 1 core(CPU)
#SBATCH --cpus-per-task=1
#SBATCH --job-name=PLINK_array   # sensible name for the job
#SBATCH --mem=10G                 # Default memory per CPU is 3GB.
#  SBATCH --constraint=avx2 # Use the smallmem-partition for jobs requiring < 10 GB RAM.

#SBATCH --array=0-4

#SBATCH -o ./Reports/PLINK_array_%j_%a.out    #Standar output message
#SBATCH -e ./Reports/PLINK_array_%j_%a.err    #Standar error message



# module purge

module load PLINK/1.9b_6.17-x86_64
module list 

names=("NE200_100cm_rep0" "NE200_100cm_rep1" "NE200_100cm_rep2" "NE100_100cm_rep1" "NE100_100cm_rep2")


# Get the name based on the array index

name="${names[SLURM_ARRAY_TASK_ID]}"

echo $name

# Make all the output directories 

mkdir  ./Dataset/$name/PLINK_out

mkdir  ./Dataset/$name/PLINK_out/Full_Dens
mkdir  ./Dataset/$name/PLINK_out/Half_Dens

mkdir  ./Dataset/$name/PLINK_out/Full_Dens/Default
mkdir  ./Dataset/$name/PLINK_out/Half_Dens/Default

mkdir  ./Dataset/$name/PLINK_out/Full_Dens/Norm
mkdir  ./Dataset/$name/PLINK_out/Half_Dens/Norm

mkdir  ./Dataset/$name/PLINK_out/Full_Dens/Norm_small
mkdir  ./Dataset/$name/PLINK_out/Half_Dens/Norm_small

# run the default scenario
plink --bfile ./Dataset/$name"/beforeQC" --homozyg --out ./Dataset/$name/PLINK_out/Full_Dens/Default/ROH_analyse
plink --bfile ./Dataset/$name"/beforeQC_50" --homozyg --out ./Dataset/$name/PLINK_out/Half_Dens/Default/ROH_analyse


# run the norm scenario
plink --bfile ./Dataset/$name"/beforeQC"  --homozyg --homozyg-snp 50 --homozyg-kb 500 --homozyg-density 100 --homozyg-gap 1000 --homozyg-het 0  --homozyg-window-snp 50 --homozyg-window-het 0 --homozyg-window-missing 1 --homozyg-window-threshold 0.05 --out  ./Dataset/$name/PLINK_out/Full_Dens/Norm/ROH_analyse
plink --bfile ./Dataset/$name"/beforeQC_50"  --homozyg --homozyg-snp 50 --homozyg-kb 500 --homozyg-density 100 --homozyg-gap 1000 --homozyg-het 0  --homozyg-window-snp 50 --homozyg-window-het 0 --homozyg-window-missing 1 --homozyg-window-threshold 0.05 --out  ./Dataset/$name/PLINK_out/Half_Dens/Norm/ROH_analyse


# run the norm small scenario 
plink --bfile ./Dataset/$name"/beforeQC" --homozyg --homozyg-snp 20 --homozyg-kb 500 --homozyg-density 100 --homozyg-gap 1000 --homozyg-het 0 --homozyg-window-snp 20 --homozyg-window-het 0 --homozyg-window-missing 1 --homozyg-window-threshold 0.05 --out ./Dataset/$name/PLINK_out/Full_Dens/Norm_small/ROH_analyse
plink --bfile ./Dataset/$name"/beforeQC_50" --homozyg --homozyg-snp 20 --homozyg-kb 500 --homozyg-density 100 --homozyg-gap 1000 --homozyg-het 0 --homozyg-window-snp 20 --homozyg-window-het 0 --homozyg-window-missing 1 --homozyg-window-threshold 0.05 --out ./Dataset/$name/PLINK_out/Half_Dens/Norm_small/ROH_analyse



