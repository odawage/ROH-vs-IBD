
# ===================
#   Script for calculating FIBD and FROH 
#   Both overall and for a specific BIN
# =====================

# Loading packages
# ======================
.libPaths("/mnt/users/odwa/R/myLib" )
library(tidyverse)
library(data.table)
library(foreach)
library(doSNOW)
library(doParallel)

# Cluster
# ===============
# Number of cores to use
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))
print(paste0("CPUs available: ", num_cores))
# Use one less than available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterEvalQ(cl, {
  .libPaths("/mnt/users/odwa/R/myLib" )
  library(tidyverse)
  library(data.table)
  library(foreach)
  library(doSNOW)
  library(doParallel)
})

# Combos
# ====================
proj.dir = "/mnt/project/SimData/Paper_1"

NEs = c("NE100","NE200")
Reps = c(1:10)
Densities = c("Full_Dens","Half_Dens")
Methods = c("Meyermans","Default","Norm", "Norm_small")

combinations_df <- expand.grid(NE = NEs , Rep = Reps)
combinations_all <- expand.grid(NE = NEs , Rep = Reps,Dens = Densities, Method = Methods)

# ===============
#   The same code grouped by bin as well 
#   Meaning that the split by bin code is added 
#==============

IBD_ind_all <- data.frame()


# Parallel loop     
IBD_ind_all_list <- foreach(i = 1:nrow(combinations_df), .packages = c('data.table', 'tidyverse' )) %dopar% {
  
  # Get the current combination
  NE <- combinations_df$NE[i]
  Rep <- combinations_df$Rep[i]
  
  data.dir = paste0(proj.dir,"/",NE,"/Rep",Rep)

  true_IBD_position <- read_csv(file.path(data.dir,"true_IBD_position.txt"),show_col_types = FALSE)
  IBD_geno_length <-  (true_IBD_position$position_bpp[nrow(true_IBD_position)] - true_IBD_position$position_bpp[1])/1000
  N_IBD_SNPs <- nrow(true_IBD_position)

  info_ADAM_roh_bin <- read_csv(file.path(data.dir,"info_ADAM_roh_bin.txt"),show_col_types = FALSE) %>% 
    mutate(KB = (stop_bpp - start_bpp)/1000) 
  
  # Make the BINS 
  bin_vec <- vector()
  # Define breakpoints and labels
  breakpoints <- c(0, 1000, 2000, 3000, 4000, 5000, 9000, Inf)
  labels <- c("0_1MB", "1_2MB", "2_3MB", "3_4MB", "4_5MB", "5_9MB", "9MB+")
  # Use cut to create the bin_vec
  bin_vec <- cut(info_ADAM_roh_bin$KB, breaks = breakpoints, labels = labels, right = FALSE)
  bin_vec_fac_2 <- factor(bin_vec, levels = c("0_1MB", "1_2MB" , "2_3MB" , "3_4MB" , "4_5MB","5_9MB","9MB+" ))
  info_ADAM_roh_bin$BIN <- bin_vec_fac_2

  # Group and summarise
  F_tot <- info_ADAM_roh_bin %>% 
    group_by(ID) %>% 
    dplyr::summarise(IBD_SNPs = sum(NSNP),
                     sum_IBD_kb = sum(KB),
                     F_IBD_snp = IBD_SNPs/N_IBD_SNPs,
                     F_IBD_kb = sum_IBD_kb/IBD_geno_length) %>% 
    mutate(BIN = "Total",
           NE = NE,
           Rep = Rep)
  
  IBD_ind <- info_ADAM_roh_bin %>% 
    group_by(ID, BIN) %>% 
    dplyr::summarise(IBD_SNPs = sum(NSNP),
              sum_IBD_kb = sum(KB),
              F_IBD_snp = IBD_SNPs/N_IBD_SNPs,
              F_IBD_kb = sum_IBD_kb/IBD_geno_length) %>% 
    mutate(NE = NE,
           Rep = Rep)
  IBD_ind <-rbind(IBD_ind, F_tot)
  return(IBD_ind)
  #cat(dataset, "\n")
} 

IBD_ind_all <- do.call(rbind, IBD_ind_all_list)

# ==============
#   Rep BIN FROH 
# ===============


ROH_ind_all_list <- foreach(i = 1:nrow(combinations_all), .packages = c('data.table', 'tidyverse' )) %dopar% {

  NE = combinations_all$NE[i]
  Rep = combinations_all$Rep[i]
  Dens = combinations_all$Dens[i]
  Method = combinations_all$Method[i]
  
  dens.dir = paste0(proj.dir,"/",NE,"/Rep",Rep,"/",Dens)
  roh.dir = paste0(proj.dir,"/",NE,"/Rep",Rep,"/",Dens,"/",Method)
  
  markerPosition_kb <- read_csv(file.path(dens.dir,"markerPositionsRep1.res.bz2"),show_col_types = FALSE) %>%
    mutate(position_kb = BPP/1000)
  geno_length <-  markerPosition_kb$position_kb[nrow(markerPosition_kb)] - markerPosition_kb$position_kb[1]
  N_SNPs <- nrow(markerPosition_kb)

  
  ROH_edit <- fread(file.path(roh.dir,"ROH_analyse.hom")) %>% 
    dplyr::rename(ID = "IID")
  
  # Make the BINS
  bin_vec <- vector()
  # Define breakpoints and labels
  breakpoints <- c(0, 1000, 2000, 3000, 4000, 5000, 9000, Inf)
  labels <- c("0_1MB", "1_2MB", "2_3MB", "3_4MB", "4_5MB", "5_9MB", "9MB+")
  # Use cut to create the bin_vec
  bin_vec <- cut(ROH_edit$KB, breaks = breakpoints, labels = labels, right = FALSE)
  bin_vec_fac_2 <- factor(bin_vec, levels = c("0_1MB", "1_2MB" , "2_3MB" , "3_4MB" , "4_5MB","5_9MB","9MB+" ))
  ROH_edit$BIN <- bin_vec_fac_2
  
  
  # summarise
  F_tot <- ROH_edit %>% 
    group_by(ID) %>% 
    dplyr::summarise(ROH_SNPs = sum(NSNP),
                     sum_ROH_kb = sum(KB),
                     F_ROH_snp = ROH_SNPs/N_SNPs,
                     F_ROH_kb = sum_ROH_kb/geno_length) %>% 
    mutate(BIN = "Total",
           NE = NE,
           Rep = Rep, 
           Dens = Dens,
           Method = Method)
  
  ROH_ind <- ROH_edit %>% 
    group_by(ID, BIN) %>% 
    dplyr::summarise(ROH_SNPs = sum(NSNP),
              sum_ROH_kb = sum(KB),
              F_ROH_snp = ROH_SNPs/N_SNPs,
              F_ROH_kb = sum_ROH_kb/geno_length) %>% 
    mutate(NE = NE,
           Rep = Rep, 
           Dens = Dens,
           Method = Method)
  
  ROH_ind <- rbind(ROH_ind,F_tot)
  return(ROH_ind) 
} #set_dens

ROH_ind_all <- do.call(rbind, ROH_ind_all_list)


# Here is the final funny little part
# So since not all IBD segments are picked up, 
# and not all ROH are IBD there is a section of animals missing combinations of 
# parametersettings, density, BIN etc
# So here I use the Complete() to make all the missing combinations 




tester <- ROH_ind_all %>% 
  group_by(NE,Rep,Dens,Method) %>% 
  complete(ID, BIN ) %>% 
  mutate_at(c(1:9), ~ replace_na(.,0))

#then join it all together 
corr_basis_BIN  <- left_join( tester,IBD_ind_all, by = c("ID", "BIN","NE","Rep")) %>% 
  replace(is.na(.), 0)

write.table(x = corr_basis_BIN ,file= file.path(proj.dir,"Results/F_ROH_IBD_BIN.txt"),row.names=FALSE,quote=FALSE)




