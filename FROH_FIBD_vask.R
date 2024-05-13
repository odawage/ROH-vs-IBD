
# ===================
#   Script for calculating FIBD and FROH 
#   Both overall and for a specific BIN
# =====================


# Loading packages
# ======================
library(tidyverse)
library(data.table)

# Loading the scenario list
# ====================

scenario_tbl <- read_table("Scenarios",
                           col_names = FALSE)[1:40,]
colnames(scenario_tbl) <- c("Dataset","Density","Setting")

datasets <- unique(scenario_tbl$Dataset)



# Calculating per individual FROH
# =================

ROH_ind_all <- data.frame()
IBD_ind_all <- data.frame()

# I'm looping over dataset cause it 
# it's the remnants of previous version 
for (dataset in unique(scenario_tbl$Dataset)  ) {

  temp_scen <- scenario_tbl %>% 
    filter(Dataset == dataset)
  
  
  for (i in 1:8) {
    setting <- temp_scen$Setting[i]
    density <- temp_scen$Density[i]
    
    # To get the right proportion of SNPs I have to choose different marker maps
    if (density == "Half") {
     markerPosition_kb <- read_csv(paste0("Dataset/",dataset,"/50_markerPositionsRep.res.bz2"),show_col_types = FALSE) %>%
       mutate(position_kb = BPP/1000)
    geno_length <-  markerPosition_kb$position_kb[nrow(markerPosition_kb)] - markerPosition_kb$position_kb[1]
    N_SNPs <- nrow(markerPosition_kb)
    }else{
      markerPosition_kb <- read_csv(paste0("Dataset/",dataset,"/R_output/markerPosition_kb.txt"),show_col_types = FALSE)%>%
        mutate(position_kb = position_bpp/1000)
      
      geno_length <-  markerPosition_kb$position_kb[nrow(markerPosition_kb)] - markerPosition_kb$position_kb[1]
      N_SNPs <- nrow(markerPosition_kb)
    }
    
    
    ROH_edit <- read_csv(paste0("Dataset/",dataset,"/R_output/",density,"_Dens/",setting,"/ROH_edit.txt"),show_col_types = FALSE)
    
    ROH_ind <- ROH_edit %>% 
      group_by(ID) %>% 
      summarise(ROH_SNPs = sum(NSNP),
                sum_kb = sum(KB),
                F_ROH_snp = ROH_SNPs/N_SNPs,       # This is the one used in the Paper 
                F_ROH_kb = sum_kb/geno_length) %>% 
      mutate(dataset = dataset, # Showing in the scenario information 
             density = density,
             setting = setting)
    
    ROH_ind_all <- rbind(ROH_ind_all, ROH_ind)
    
  } #set_dens
  cat(dataset, "\n")
}



write.table(x=ROH_ind_all ,file= paste("./Dataset/Correlation/F_ROH.txt", sep = ""),row.names=FALSE,quote=FALSE)


# ==========
#
# Then repeat for FIBD
#
# =============

IBD_ind_all <- data.frame()

for (dataset in unique(scenario_tbl$Dataset)  ) {
  
  
  true_IBD_position <- read_csv(paste0("Dataset/",dataset,"/R_output/true_IBD_position.txt"),show_col_types = FALSE)
  IBD_geno_length <-  true_IBD_position$position_bpp[nrow(true_IBD_position)] - true_IBD_position$position_bpp[1]
  N_IBD_SNPs <- nrow(true_IBD_position)
  
  
  info_ADAM_roh_bin <- read_csv(paste0("Dataset/",dataset,"/IBD_segments/info_ADAM_roh_bin.txt"),show_col_types = FALSE) %>% 
    mutate(NSNPS = stop -start, 
           BPP = stop_bpp - start_bpp) 
  
  IBD_ind <- info_ADAM_roh_bin %>% 
    group_by(ID) %>% 
    summarise(IBD_SNPs = sum(NSNPS),
              sum_IBD_BPP = sum(BPP),
              F_IBD_snp = IBD_SNPs/N_IBD_SNPs, # Useing this one in the paper
              F_IBD_BPP = sum_IBD_BPP/IBD_geno_length) %>% 
    mutate(dataset = dataset)
  
  IBD_ind_all <- rbind(IBD_ind_all, IBD_ind)
}


write.table(x= IBD_ind_all ,file= paste("./Dataset/Correlation/F_IBD.txt", sep = ""),row.names=FALSE,quote=FALSE)




# ===============
#   The same code grouped by bin as well 
#   Meaning that the split by bin code is added 
#==============



ROH_ind_all <- data.frame()
IBD_ind_all <- data.frame()

for (dataset in unique(scenario_tbl$Dataset)  ) {
  
  
  true_IBD_position <- read_csv(paste0("Dataset/",dataset,"/R_output/true_IBD_position.txt"),show_col_types = FALSE)
  IBD_geno_length <-  (true_IBD_position$position_kb[nrow(true_IBD_position)] - true_IBD_position$position_kb[1])/1000
  N_IBD_SNPs <- nrow(true_IBD_position)
  
  
  
  info_ADAM_roh_bin <- read_csv(paste0("Dataset/",dataset,"/IBD_segments/info_ADAM_roh_bin.txt"),show_col_types = FALSE) %>% 
    mutate(NSNPS = stop -start, 
           KB = (stop_bpp - start_bpp)/1000) 
  
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
  IBD_ind <- info_ADAM_roh_bin %>% 
    group_by(ID, BIN) %>% 
    summarise(IBD_SNPs = sum(NSNPS),
              sum_IBD_kb = sum(KB),
              F_IBD_snp = IBD_SNPs/N_IBD_SNPs,
              F_IBD_kb = sum_IBD_kb/IBD_geno_length) %>% 
    mutate(dataset = dataset)
  
  IBD_ind_all <- rbind(IBD_ind_all, IBD_ind)
  
  
  temp_scen <- scenario_tbl %>% 
    filter(Dataset == dataset)
  
  
  for (i in 1:8) {
    setting <- temp_scen$Setting[i]
    density <- temp_scen$Density[i]
    
    if (density == "Half") {
      markerPosition_kb <- read_csv(paste0("Dataset/",dataset,"/50_markerPositionsRep.res.bz2"),show_col_types = FALSE) %>%
        mutate(position_kb = BPP/1000)
      geno_length <-  markerPosition_kb$position_kb[nrow(markerPosition_kb)] - markerPosition_kb$position_kb[1]
      N_SNPs <- nrow(markerPosition_kb)
    }else{
      markerPosition_kb <- read_csv(paste0("Dataset/",dataset,"/R_output/markerPosition_kb.txt"),show_col_types = FALSE)%>%
        mutate(position_kb = position_bpp/1000)
      
      geno_length <-  markerPosition_kb$position_kb[nrow(markerPosition_kb)] - markerPosition_kb$position_kb[1]
      N_SNPs <- nrow(markerPosition_kb)
    }
    
    
    ROH_edit <- read_csv(paste0("Dataset/",dataset,"/R_output/",density,"_Dens/",setting,"/ROH_edit.txt"),show_col_types = FALSE)
    
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
    ROH_ind <- ROH_edit %>% 
      group_by(ID, BIN) %>% 
      summarise(ROH_SNPs = sum(NSNP),
                sum_kb = sum(KB),
                F_ROH_snp = ROH_SNPs/N_SNPs,
                F_ROH_kb = sum_kb/geno_length) %>% 
      mutate(dataset = dataset, 
             density = density,
             setting = setting)
    
    ROH_ind_all <- rbind(ROH_ind_all, ROH_ind)
    
  } #set_dens
  cat(dataset, "\n")
}


# Here is the final funny little part
# So since not all IBD segments are picked up, 
# and not all ROH are IBD there is a section of animals missing combinations of 
# parametersettings, density, BIN etc
# So here I use the Complete() to make all the missing combinations 

tester <- ROH_ind_all %>% 
  group_by(dataset,density,setting ) %>% 
  complete(ID, BIN ) %>% 
  mutate_at(c(1:9), ~replace_na(.,0))

#then join it all together 
corr_basis_BIN  <- left_join( tester,IBD_ind_all, by = c("ID", "BIN","dataset")) %>% 
  replace(is.na(.), 0)

write.table(x = corr_basis_BIN ,file= paste("./Dataset/Correlation/F_ROH_IBD_BIN.txt", sep = ""),row.names=FALSE,quote=FALSE)




