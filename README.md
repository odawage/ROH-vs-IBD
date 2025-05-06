# ROH-vs-IBD
Code related to paper ""

# General info 
  The scripts are written for SLURM, R 4.3.3 and PLINK 1.9 and 2.0
  The "Scenarios" file is a list of all paramenter combinations, dataset, and marker densities. 

# Folder structure for data
Effective population size
 └── Data of replicat "X"
     └── Density 

The rest of the result folders are created by the scripts

# Order of running the files :yellow_heart:

1. Raw_to_PedMap.R
   - Convertes files from ADAMs format to .ped/.map and BIM, BED, FAM
     It also makes the split between test(masked) and train markers
     You have to choose the relevant generations 

3. Detect_ROH.R
   - Runs the script with the ROH detection scenarios from v2_ROH_detection_functions.R

5. IBD_segments.R
   - finds the IBD segments per dataset.
     Also select relevant generations here.
  
7.  ROH_vs_IBD.R
   - Calculated rate of true positive and power for all scenarios

9.  R_F_in_ROH.R
    - All steps needed to find F|ROH for all scenarios 

11. FROH_FIBD.R
    - Finds FROH and FIBD for all animals.
      One overall, and one per BIN estimate.

# Licensing 

Creative commons. 🫶
