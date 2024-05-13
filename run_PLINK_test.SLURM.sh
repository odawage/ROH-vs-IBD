#!/bin/bash
#SBATCH --ntasks=1              # 1 core(CPU)
#SBATCH --cpus-per-task=1
#SBATCH --job-name=PLINK_array   # sensible name for the job
#SBATCH --mem=3G                 # Default memory per CPU is 3GB.
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

mkdir  ./Dataset/$name/PLINK_out

mkdir  ./Dataset/$name/PLINK_out/Full_Dens
mkdir  ./Dataset/$name/PLINK_out/Half_Dens

mkdir  ./Dataset/$name/PLINK_out/Full_Dens/Tester
mkdir  ./Dataset/$name/PLINK_out/Half_Dens/Tester

#plink --bfile ./Dataset/$name"/beforeQC" --homozyg --homozyg-window-snp 100 --out ./Dataset/$name#/PLINK_out/Full_Dens/Tester/ROH_analyse
#plink --bfile ./Dataset/$name"/beforeQC_50" --homozyg --homozyg-window-snp 100 --out ./Dataset/$name/PLINK_out/Half_Dens/Tester/ROH_analyse

plink --bfile ./Dataset/$name"/beforeQC" --homozyg  ---homozyg-snp 50 --homozyg-window-snp 50 --homozyg-kb 500 --homozyg-window-threshold 0.1  --out ./Dataset/$name/PLINK_out/Full_Dens/Tester/ROH_analyse
plink --bfile ./Dataset/$name"/beforeQC_50" --homozyg  --homozyg-snp 50 --homozyg-window-snp 50 --homozyg-kb 500 --homozyg-window-threshold 0.1 --out ./Dataset/$name/PLINK_out/Half_Dens/Tester/ROH_analyse
