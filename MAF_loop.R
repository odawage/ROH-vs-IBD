

# 
# source("/usr/share/lmod/lmod/init/R")
# Sys.setenv(MODULEPATH = '/cluster/modules/all')
# module("load PLINK/1.9b_6.17-x86_64")
# plink = "plink"


# =======================
# Libraries 
.libPaths("/mnt/users/odwa/R/myLib" )


library(data.table)
library(plyr)
library(dplyr)
#install.packages('openxlsx')
library(openxlsx)
library(openxlsx)
library(ggplot2)
library(readr)
library(stringi)

library(doParallel)
library(foreach)

source("~/paper-1/ROH-vs-IBD/containment_zone/ROH_detection_functions.R")
source("~/paper-1/ROH-vs-IBD/containment_zone/vcferr_spTallman.R")
source("~/paper-1/ROH-vs-IBD/containment_zone/IBD_vs_ROH_functions.R")



input.dir <- "/mnt/users/odwa/PLINK/Code/Dataset"
#dataset = "NE100_100cm_rep1"
inped = "beforeQC"
tmp.dir = "/net/fs-2/scale/OrionStore/Scratch/odwa/paper-1"
datasets = c("NE200_100cm_rep0", "NE200_100cm_rep1", "NE200_100cm_rep2", "NE100_100cm_rep1", "NE100_100cm_rep2")


# dataset = "NE100_100cm_rep1"
# inped = "beforeQC"
# dens = ""
# MAF = 0

# 
# 
# for (dataset in datasets) {
#   for (dens in c("","_50")){
#     for (MAF in c(0 , 0.01,0.05)) {
#       if (grepl("_50", dens)) {
#         dens.dir <-paste0(tmp.dir,"/",dataset,"/Half_Dens") 
#       } else {
#         dens.dir <-paste0(tmp.dir,"/",dataset,"/Full_Dens") 
#       }
#       
#       dir.create(dens.dir)
#     
#       if (MAF == 0) {
#       system(paste0("plink --bfile ",input.dir,"/",dataset,"/",inped,dens,
#                       " --make-bed --out ",dens.dir,"/inped_maf"))
#       }else{
#       system(paste0("plink --bfile ",input.dir,"/",dataset,"/",inped,dens,
#                     " --maf ",MAF," --make-bed --out ",dens.dir,"/inped_maf"))
#       } 
#     
#       # Detect the ROH
#       Gorssen(input.dir = dens.dir, dataset = dataset, paste0("inped_maf"), tmp.dir = tmp.dir, dens = dens, err = paste0(MAF,"_MAF"))
#       Default(input.dir = dens.dir, dataset = dataset, paste0("inped_maf"), tmp.dir = tmp.dir, dens = dens, err = paste0(MAF,"_MAF"))
#       Norm(input.dir = dens.dir, dataset = dataset, paste0("inped_maf"), tmp.dir = tmp.dir, dens = dens, err = paste0(MAF,"_MAF"))
#       Norm_small(input.dir = dens.dir, dataset = dataset, paste0("inped_maf"), tmp.dir = tmp.dir, dens = dens, err = paste0(MAF,"_MAF"))
#       
#     } # err
#   } # dens
# } # dataset
# 
# 
# for (dataset in datasets) {
#   for (dens in c("","_50")){
#     for (MAF in c(0 , 0.01,0.05)) {
#       for (gen in c(20,40,60,80,100)) {
#         
#       
#       if (grepl("_50", dens)) {
#         dens.dir <-paste0(tmp.dir,"/",dataset,"/Half_Dens") 
#       } else {
#         dens.dir <-paste0(tmp.dir,"/",dataset,"/Full_Dens") 
#       }
#       
#       dir.create(dens.dir)
#       
#       fread(paste0(input.dir,"/",dataset,"/pedigree_rel.generations.txt")) %>% 
#         filter(birth == gen) %>% 
#         mutate(FID = 0) %>% 
#         select(FID,id) %>% 
#         write.table(file = paste0(dens.dir,"/rel_inds.txt"), 
#                     row.names = F, 
#                     col.names = F, 
#                     quote = F, 
#                     sep = '\t')
#       
#       if (MAF == 0) {
#           
#         system(paste0("plink --bfile ",input.dir,"/",dataset,"/",inped,dens,
#                       " --keep ",dens.dir,"/rel_inds.txt --make-bed --out ",dens.dir,"/inped_maf"))
#       }else{
#         system(paste0("plink --bfile ",input.dir,"/",dataset,"/",inped,dens,
#                       " --keep ",dens.dir,"/rel_inds.txt --maf ",MAF," --make-bed --out ",dens.dir,"/inped_maf"))
#       } 
#       
#       # Detect the ROH
#       Gorssen(input.dir = dens.dir, dataset = dataset, paste0("inped_maf"), tmp.dir = tmp.dir, dens = dens, err = paste0(MAF,"_",gen,"_MAF"))
#       Default(input.dir = dens.dir, dataset = dataset, paste0("inped_maf"), tmp.dir = tmp.dir, dens = dens, err = paste0(MAF,"_",gen,"_MAF"))
#       Norm(input.dir = dens.dir, dataset = dataset, paste0("inped_maf"), tmp.dir = tmp.dir, dens = dens, err = paste0(MAF,"_",gen,"_MAF"))
#       Norm_small(input.dir = dens.dir, dataset = dataset, paste0("inped_maf"), tmp.dir = tmp.dir, dens = dens, err = paste0(MAF,"_",gen,"_MAF"))
#       } # generation 
#     } # err
#   } # dens
# } # dataset
# 






# ===================
#
#   Benchmarking 
#
# ===================

# Number of cores to use
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

# Use one less than available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the different values
densities <- c("Full_Dens", "Half_Dens")
methods <- c("Meyermans", "Default", "Norm", "Norm_small")
MAF <- c(0, 0.01,0.05)

# Generate all combinations of dataset, density, method, and err
combinations <- expand.grid(dataset = datasets, dens = densities, method = methods, MAF = MAF)

# Parallel processing over all combinations

# Load required packages on each worker
clusterEvalQ(cl, {
  .libPaths("/mnt/users/odwa/R/myLib" )
  library(data.table)
  library(plyr)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  library(readr)
  library(stringi)
  library(foreach)
  library(tidyr)
})

# Load custom functions for each worker
clusterEvalQ(cl, {
  source("~/paper-1/ROH-vs-IBD/containment_zone/ROH_detection_functions.R")
  source("~/paper-1/ROH-vs-IBD/containment_zone/vcferr_spTallman.R")
  source("~/paper-1/ROH-vs-IBD/containment_zone/IBD_vs_ROH_functions.R")
})

# Parallel loop
foreach(i = 1:nrow(combinations),.packages = c('data.table', 'plyr', 'dplyr', 'readr',"tidyr")) %dopar% {
  
  # Get the current combination
  dataset <- combinations$dataset[i]
  dens <- combinations$dens[i]
  method <- combinations$method[i]
  MAF <- combinations$MAF[i]
  
  data.input.dir <- paste0(input.dir, "/", dataset)
  
  marker_pos <- read_csv(paste0(data.input.dir, "/R_output/markerPosition_kb.txt"), show_col_types = FALSE)
  true_IBD <- read_csv(paste0(data.input.dir, "/R_output/All_Ind_true_IBD.res.bz2"), show_col_types = FALSE)
  IBD_pos <- read.csv(paste0(data.input.dir, "/R_output/true_IBD_position.txt"))
  ped <- read_table(paste0(data.input.dir, "/pedigree_rel.generations.txt"))
  obs_IBD <- read_csv(paste(data.input.dir, "/IBD_segments/extend_bin_ADAM_ROH.res.bz2", sep = ""))
  obs_homo <- read_csv(paste0(data.input.dir, "/R_output/all_ind_observed_homozyg_QTL.res.bz2"))
  
  info_ADAM_ROH_bin <- read.csv(paste0(data.input.dir, "/IBD_segments/info_ADAM_roh_bin.txt")) %>%
    filter(ID %in% ped$id) %>%
    mutate(start_kb = start_bpp / 1000,
           stop_kb = stop_bpp / 1000)
  
  ROH_inf <- data.frame()
  for (gen in c(20,40,60,80,100)) {
    roh.dir <- paste0(tmp.dir, "/", dataset, "/", dens, "/", method, "/", MAF,"_",gen,"_MAF")
    
    ROH_inf_t <- read_table(paste0(roh.dir, "/ROH_analyse.hom")) %>%
      select(c(IID, SNP1, SNP2, POS1, POS2, KB, NSNP)) %>%
      rename(ID = IID) %>%
      filter(ID %in% ped$id) %>%
      group_by(ID) %>%
      mutate(run = row_number())
    if (nrow(ROH_inf_t) > 0) {
      ROH_inf <- rbind(ROH_inf, ROH_inf_t)  
    }
    
  }
  
  PLINK_ROH_bin <- RohBinGeno(ROH_inf, marker_pos)
  
  RohPower(info_ADAM_ROH_bin, PLINK_ROH_bin, marker_pos, tmp.dir = roh.dir)
  RohTruePos(ROH_inf, obs_IBD, obs_homo, marker_pos, tmp.dir = roh.dir)
  
} # End of foreach loop

# Stop the cluster after processing
stopCluster(cl)
