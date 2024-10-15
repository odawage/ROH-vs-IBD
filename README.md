# ROH-vs-IBD
Code related to paper ""

# General info 
  The scripts are written for SLURM, R 4.3.3 and PLINK 1.9 and 2.0
  The "Scenarios" file is a list of all paramenter combinations, dataset, and marker densities. 

# Folder structure 

- Code
  - Dataset
    - [name of dataset]
      - [set all raw ADAM files here]

The rest of the result folders are created by the scripts

# Order of running the files :yellow_heart:

1. Raw_to_PedMap.R
   - Convertes files from ADAMs format to .ped/.map and BIM, BED, FAM
     It also makes the split between test(masked) and train markers
     You have to choose the relevant generations 

3. run_Gorrsen_V3.SLURM.sh and Gorrsen_V3.R
   - Runs the script for finding ROH using Meyermans formulas. 

3. run_PLINK_parameters.SLURM.sh
   - Running parameters for default, norm, and norm small in PLINK.
     runs for both full and half density markers 

5. Run_IBD_segments.sh & IBD_segments.R
   - finds the IBD segments per dataset.
     Also select relevant generations here.
  
7. run_ROH_vs_IBD.sh  & ROH_vs_IBD.R
   - Calculated rate of true positive and power for all scenarios

9. run_F_in_ROH & R_F_in_ROH.R
    - All steps needed to find F|ROH for all scenarios 

11. FROH_FIBD.R
    - Finds FROH and FIBD for all animals.
      One overall, and one per BIN estimate.

# Licensing 

Creative commons. Live, laugh, love 🫶
