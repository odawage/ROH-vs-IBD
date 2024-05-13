# ROH-vs-IBD
Code related to paper ""

# General info 
  The scripts are written for SLURM, R 4.3.3 and PLINK 1.9 and 2.0
  

# ===========================
Order of running the files <3
# ==========================

1. Raw_to_PedMap.R 

  Here the files are converted from the ADAM format to .ped/.map and BIM, BED, FAM 
  It also makes the split between test and train markers
  You also have to choose the relevant generations 

2. run_array_rewrite_Gorrsen_V2.SLURM.sh and array_rewrite_Gorrsen_V2.R

  Runs the script for finding the Meyerman ROH 

3. run_PLINK_parameters.SLURM.sh 
  Running defauls, norm, and norm small 
  runs for both full and half density markers 

4. Run_IBD_segments.sh & IBD_segments.R
  finds the segments in the dataset that are the same for all settings and density 
  Also selects relevant generations here 
  
5. run_ROH_vs_IBD.sh  & ROH_vs_IBD.R 

6. run_hetero_dataset_gen & R_hetero_dataset_gen.R

7. Corr_Froh_FIBD.R
