
#libraries
.libPaths("/mnt/users/odwa/R/myLib" )
library(stringi)
library(tidyverse)
library(data.table)
library(foreach)
#library(snow)
library(doSNOW)
library(itertools)
library(doParallel)
library(parallel)
# =============
#   Directories
# =============

data.dir <- paste0("/mnt/project/SimData/Paper_1")

# =============
#  Combos 
# =============

NE = c("NE100","NE200")
Rep = c(1:10)
Dens = c("Full_Dens","Half_Dens")
Method = c("Meyermans","Default","Norm", "Norm_small")

combinations <- expand.grid(NE = NE , Rep = Rep, Dens = Dens ,Method = Method)

# ==================
#
#   Benchmarking 
#
# ===================

# Number of cores to use
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
print(paste0("CPUs available: ", num_cores))
# Use one less than available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  .libPaths("/mnt/users/odwa/R/myLib" )
  library(stringi)
  library(tidyverse)
  library(data.table)
  library(foreach)
  #library(snow)
  library(doSNOW)
  library(doParallel)
  library(parallel)
})

source("/mnt/users/odwa/paper-1/ROH-vs-IBD/IBD_vs_ROH_functions.R")

# Parallel loop
foreach(i = 1:nrow(combinations), .packages = c('data.table', 'dplyr', 'readr',"tidyr"),
        .export = c("RohBinGeno","RohPower","RohTruePos")) %dopar% {
  
  # Get the current combination
  NE <- combinations$NE[i]
  Rep <- combinations$Rep[i]
  Dens <- combinations$Dens[i]
  Method <- combinations$Method[i]

  
  data.input.dir <- paste0(data.dir, "/",NE, "/Rep",Rep)
  roh.dir <- paste0(data.input.dir,"/",Dens,"/",Method)
  
  marker_pos <- read_csv(paste0(data.input.dir, "/markerPosition_kb.txt"), show_col_types = FALSE)
  true_IBD <- read_csv(paste0(data.input.dir, "/All_Ind_true_IBD.res.bz2"), show_col_types = FALSE)
  IBD_pos <- read.csv(paste0(data.input.dir, "/true_IBD_position.txt"))
  ped <- read_table(paste0(data.input.dir, "/pedigree_rel_generations.txt"))
  obs_IBD <- read_csv(paste(data.input.dir, "/extend_bin_ADAM_ROH.res.bz2", sep = ""))
  obs_homo <- read_csv(paste0(data.input.dir, "/all_ind_observed_homozyg_QTL.res.bz2"))
  
  info_ADAM_ROH_bin <- read.csv(paste0(data.input.dir, "/info_ADAM_roh_bin.txt")) %>%
    filter(ID %in% ped$id) %>%
    mutate(start_kb = start_bpp / 1000,
           stop_kb = stop_bpp / 1000)
  

  ROH_inf <- read_table(paste0(roh.dir, "/ROH_analyse.hom")) %>%
      select(c(IID, SNP1, SNP2, POS1, POS2, KB, NSNP)) %>%
      dplyr::rename(ID = IID) %>%
      filter(ID %in% ped$id) %>%
      group_by(ID) %>%
      mutate(run = row_number())

  PLINK_ROH_bin <- RohBinGeno(ROH_inf, marker_pos)
  RohPower(info_ADAM_ROH_bin, PLINK_ROH_bin, marker_pos, tmp.dir = roh.dir)
  RohTruePos(ROH_inf, obs_IBD, obs_homo, marker_pos, tmp.dir = roh.dir)
  
} # End of foreach loop

# Stop the cluster after processing
stopCluster(cl)
