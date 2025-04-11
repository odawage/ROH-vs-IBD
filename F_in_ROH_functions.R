
####
#
#    Then All the functions needed for the looping 
#
####
process_SNP <- function(SNP, marker_pos, marker_pos_test, ROH_edit_gen, homo_bin_gen, generation) {
  marker <- marker_pos_test$locus[SNP]
  # output dataframe
  result <- data.frame(Marker = marker_pos_test$locus[SNP],
                       tmp_marker = marker_pos$tmp_marker[marker_pos$locus %in% marker],
                       Generation = generation,
                       ROH_0_1 = 0,
                       ROH_1_2 = 0,
                       ROH_2_3 = 0,
                       ROH_3_4 = 0,
                       ROH_4_5 = 0,
                       ROH_5_9 = 0,
                       ROH_9 = 0,
                       hom_ROH_0_1 = 0,
                       hom_ROH_1_2 = 0,
                       hom_ROH_2_3 = 0,
                       hom_ROH_3_4 = 0,
                       hom_ROH_4_5 = 0,
                       hom_ROH_5_9 = 0,
                       hom_ROH_9 = 0)
  # Then find all the ROH that includes the test marker 
  for (ROH in 1:nrow(ROH_edit_gen)) {
    ID <- ROH_edit_gen$IID[ROH]
    
    # If the marker is in a ROH find the homozygous binary genome for the relevant individual
    # and register if the marker is homozygous or heterozygous 
    if (ROH_edit_gen$SNP1[ROH] <= marker && marker <= ROH_edit_gen$SNP2[ROH]){
      obs_hom <- strsplit(homo_bin_gen$genome[homo_bin_gen$ID %in% ID], " ")[[1]]
      pos <- marker_pos$tmp_marker[marker_pos$locus %in% marker]
      if(obs_hom[pos] == "1" ) {  # Check for character "1"
        bin <- as.numeric(ROH_edit_gen$BIN[ROH])
        result[1, 3 + bin] <- result[1, 3 + bin] + 1
        result[1, 10 + bin] <- result[1, 10 + bin] + 1
      } else {
        bin <- as.numeric(ROH_edit_gen$BIN[ROH])
        result[1, 3 + bin] <- result[1, 3 + bin] + 1
      }
    } else {
      # Handle the else case here if needed
    }
    
  }
  return(result)
}



# so then i need a IBS in roh function 
hom_in_ROH <- function(input.dir, roh.dir, Dens, gen ) {
  
  pedigree <- read_table(paste0(input.dir,"/pedigree_rel_generations.txt"),
                         col_types = cols(sire = col_character(), dam = col_character(),
                                          X11 = col_skip()), show_col_types = FALSE)
  
  # Load SNP map 
  marker_pos <- read_csv( paste0(input.dir,"/markerPosition_kb.txt"),show_col_types = FALSE)
  # Load list of test/masked markers
  marker_pos_test <- read_table(paste0(input.dir,"/",Dens,"/test.map"),
                                col_names = FALSE,show_col_types = FALSE)
  colnames(marker_pos_test) <- c("Chr","locus", "cM", "BPP")
  
  # Load homozygous binary genome
  homo_bin <- read_csv(paste0(input.dir,"/all_ind_observed_homozyg_QTL.res.bz2"),show_col_types = FALSE)
  # Load list of ROH
  ROH_edit <- read_table(paste0(roh.dir,"/ROH_analyse.hom"),show_col_types = FALSE)
  
  # Categorise the ROH based on what BIN they are in
  bin_vec <- vector()
  breakpoints <- c(0, 1000, 2000, 3000, 4000, 5000, 9000, Inf)
  labels <- c("0_1MB", "1_2MB", "2_3MB", "3_4MB", "4_5MB", "5_9MB", "9MB+")
  # Use cut to create the bin_vec
  bin_vec <- cut(ROH_edit$KB, breaks = breakpoints, labels = labels, right = FALSE)
  bin_vec <- factor(bin_vec, levels = c("0_1MB", "1_2MB" , "2_3MB" , "3_4MB" , "4_5MB","5_9MB","9MB+" ))
  ROH_edit$BIN <- bin_vec
  
  # Initialize all_result with the desired number of columns
  all_result <- data.frame( Marker = numeric(0),
                            tmp_marker= numeric(0),
                            Generation = numeric(0),
                            ROH_0_1 = numeric(0),
                            ROH_1_2 = numeric(0),
                            ROH_2_3 = numeric(0),
                            ROH_3_4 = numeric(0),
                            ROH_4_5 = numeric(0),
                            ROH_5_9 = numeric(0),
                            ROH_9 = numeric(0),
                            hom_ROH_0_1 = numeric(0),
                            hom_ROH_1_2 = numeric(0),
                            hom_ROH_2_3 = numeric(0),
                            hom_ROH_3_4 = numeric(0),
                            hom_ROH_4_5 = numeric(0),
                            hom_ROH_5_9 = numeric(0),
                            hom_ROH_9 = numeric(0))
  
  
  gen_ID <- pedigree %>%
    filter(birth %in% gen) %>%
    select(id)
  
  homo_bin_gen <- homo_bin[homo_bin$ID %in% gen_ID$id,]
  ROH_edit_gen <- ROH_edit[ROH_edit$IID %in% gen_ID$id, ]
  
  # List of SNPs to process
  SNPs <- seq_along(marker_pos_test$locus)
  
  # Use foreach to parallelize the loop over SNPs
  all_result = data.frame()
  for (SNP in SNPs) {
    result <- process_SNP(SNP, marker_pos, marker_pos_test, ROH_edit_gen, homo_bin_gen, gen)
    all_result = rbind(all_result,result) 
  }
  
  write.table(all_result, paste0(roh.dir,"/",gen,"_MAF_ROH_state_F_ROH_bin.txt"),
              row.names = FALSE, quote = FALSE)
} # scen


#####
#
#
#
####

process_ID <- function(i, marker_pos_test, ROH_edit_gen, homo_bin_gen) {
  obs_hom <- strsplit(homo_bin$genome[i], " ")[[1]]
  
  result <- data.frame(locus = marker_pos_test$locus,
                       tmp_marker = marker_pos_test$tmp_marker,
                       Generation = homo_bin$Generation[i],
                       ID = homo_bin$ID[i],
                       obs_hom = 0)
  
  vec_hom <- numeric(0)
  for (marker in marker_pos_test$tmp_marker) {
    if (obs_hom[marker] == "1") {
      vec_hom <- c(vec_hom, 1)
    }else{
      vec_hom <- c(vec_hom, 0)
    }# else
    
    
  } # test_marker
  result$obs_hom <- vec_hom
  return(result)
}


F_in_ROH <- function(homozyg_marker_all, marker_freq_all, input.dir, roh.dir,Gen){
  

  ROH_state <- read_table(paste0(roh.dir,"/",Gen,"_MAF_ROH_state_F_ROH_bin.txt"),show_col_types = FALSE)
  # Rename from so a string pattern can be used to pivot longer
  colnames(ROH_state) <- c("Marker","tmp_marker", "Generation",
                           "ROH-0_1", "ROH-1_2", "ROH-2_3" ,"ROH-3_4", "ROH-4_5", "ROH-5_9", "ROH-9",
                           "hom_ROH-0_1", "hom_ROH-1_2", "hom_ROH-2_3", "hom_ROH-3_4", "hom_ROH-4_5",
                           "hom_ROH-5_9", "hom_ROH-9")
  
  ROH_state <- ROH_state %>%
    rowwise()  %>%
    dplyr::mutate(`ROH-tot` = sum(c_across(cols = starts_with('ROH'))),
           `hom_ROH-tot` = sum(c_across(cols = starts_with('hom_ROH'))) )
  
  ROH_state <- ROH_state %>%
    pivot_longer(
      cols = starts_with("ROH") | starts_with("hom_ROH"),
      names_to = c(".value", "BIN"),
      names_sep = "-"
    ) %>%
    dplyr::rename(
      hom_obs = hom_ROH,
      ROH_obs = ROH
    )
  # attach frequency information for the test markers
  ROH_state_freq <- left_join(ROH_state, marker_freq_all[,c(2,5,7)], by = c("Marker" = "SNP"),relationship = "many-to-many")
  # attach homozyg count 
  ROH_state_freq_2 <- left_join(ROH_state_freq, homozyg_marker_all[, c(2,4:8)],
                                by = c("Marker" = "locus"),
                                relationship = "many-to-many")
  
  # Since the masked markers were selected based on MAF in generation 100, 
  # I do a second MAF filtering here to filter all crazy markers
  tester = ROH_state_freq_2 %>%
    filter(ROH_obs > 0,
           MAF > 0.05) %>%
    group_by(Generation, BIN) %>%
    dplyr::summarise(sum_ROH = sum(ROH_obs),
                     sum_hom = sum(hom_obs),
                     sum_hom_ROH = sum_hom/sum_ROH,
                     mean_2pq = mean(2*MAF* (1-MAF)),
                     mean_hom_frac = mean(obs_homo/n_ind)
    )
  
  # F_ROH_sum_frq = based on marker freq. not used in the paper
  # F_ROH_sum_frac = based on fraction of homosyg. Used in the paper
  t2 = tester %>%
    mutate(F_ROH_sum_frq = 1 -( (1-sum_hom_ROH )/ mean_2pq),
           F_ROH_sum_frac = 1 -( (1-sum_hom_ROH )/ mean_hom_frac)
    )
  return(t2)
} #scenario



# ============
#
#
# ==============




Allele_freq_gen <- function(input.path, rel_gen) {
  
  full_ped <- read_csv(file.path(input.path,"/full_split.ped"),
                       col_names = FALSE,show_col_types = FALSE)
  colnames(full_ped) <- c("ID", "Geno")
  
  
  pedigree <- read_table(file.path(input.path,"/pedigree_rel_generations.txt"),
                         col_types = cols(birthHerd = col_skip(),
                                          polyTbv1 = col_skip(), qtlTbv1 = col_skip(),
                                          residual1 = col_skip(), X11 = col_skip()),show_col_types = FALSE)
  
  dir.create(file.path(input.path,"AlleleFreq"))
  dir.create(file.path(input.path,"AlleleFreq/temp"))
  temp.dir <- file.path(input.path,"AlleleFreq/temp")
  
  tot_AlFreq <- data.frame("#CHROM" = numeric(0), "ID" = numeric(0), "REF" = numeric(0), "ALT" = numeric(0) , "ALT_FREQS" = numeric(0),"OBS_CT" = numeric(0), GEN = numeric(0))
  # then calculate the generation relative allele freq
  for (i in rel_gen) {
    gen <- paste(i)
    gen_ID <- pedigree %>%
      filter(birth %in% i) %>%
      select(id)
    gen_i_ped <- full_ped %>%
      filter(ID %in% gen_ID$id)
    
    # Make the .ped per generation 
    info <- data.frame(matrix(nrow = nrow(gen_i_ped), ncol = 6, 0))
    info[ , 2] <- as.numeric(gen_i_ped$ID)
    info <- unite(info, info, c(1:6), sep = " ")
    
    info$Geno <- gen_i_ped$Geno
    write.table(info, file= paste0(temp.dir,"/",i,"_full_gen.ped"),row.names = FALSE, col.names = FALSE,sep = " ",quote = FALSE )
    
    # make binary fileset
    system(paste0("plink --ped ",temp.dir,"/",i,"_full_gen.ped",
                  "  --map ",input.path,"/Full_Dens/full.map",
                  "  --make-bed --out ",temp.dir,"/",i,"_full_gen"))
    
    # Then we calc the freqs for all the generations
    system(paste0("plink --bfile ",temp.dir,"/",i,"_full_gen",
                  " --freq --out ",temp.dir,"/",i,"_full_gen"))
    
    # Read and combine it all into one dataframe
    Allele_freq <- read_table(paste0(temp.dir,"/",i,"_full_gen.frq") )
    Allele_freq$GEN <- seq(i,i, length.out = nrow(Allele_freq) )
    
    tot_AlFreq <- rbind(tot_AlFreq, Allele_freq)
  } # rel gen
  
  tot_AlFreq_data <- tot_AlFreq %>%
    mutate(NE = NE, 
           Rep = Rep)
  
  print(unlink(temp.dir, recursive=TRUE))
  
  # Write that dataframe
  write.table(tot_AlFreq_data, file= file.path(input.path,"AlleleFreqs_gens.txt"),sep = "," , quote = FALSE, row.names = FALSE)
}




# ==================
#
# Now count the hom per gen per dataset regardless of masked marker being in a ROH or not
#
# ==================


# ==================
# The function we're feeding to the cores
# This one is running parallell over individual 
# ==================

process_homo <- function(i, marker_pos_test, homo_bin) {
  obs_hom <- strsplit(homo_bin$genome[i], " ")[[1]]
  
  result <- data.frame(locus = marker_pos_test$locus,
                       tmp_marker = marker_pos_test$tmp_marker,
                       Generation = homo_bin$Generation[i],
                       ID = homo_bin$ID[i],
                       obs_hom = 0)
  
  vec_hom <- numeric(0)
  for (marker in marker_pos_test$tmp_marker) {
    if (obs_hom[marker] == "1") {
      vec_hom <- c(vec_hom, 1)
    }else{
      vec_hom <- c(vec_hom, 0)
    }# else
    
    
  } # test_marker
  result$obs_hom <- vec_hom
  return(result)
}


# ===========
# The nested loop over dataset and density 
# The different densities have different masked markers
# even though they are in the same dataset. 
# ===========


homo_test <- function(data.dir, NE, Rep, Dens, rel_gen) {
  # The only thing I need here is what year the individual is born in
  pedigree <- read_table(file.path(data.dir,"pedigree_rel_generations.txt"),
                         col_types = cols(sire = col_character(), dam = col_character(),
                                          X11 = col_skip()), show_col_types = FALSE)
  
  # The pre-processed homo/hetero state for all the markers
  homo_bin <- read_csv(file.path(data.dir,"all_ind_observed_homozyg_QTL.res.bz2"), show_col_types = FALSE)
  # attach the generation, so I only need to search one dataframe in the ID loop
  homo_bin <- left_join(homo_bin, pedigree[,c(1,5)], by = c("ID" = "id") ) %>%
    mutate(Generation = birth)
  
  # since I had to remove markers that didn't have any position i needed to create a second
  # set of names. So TMP_marker is the observable position in the genome.
  marker_pos <- read_csv(file.path(data.dir,"markerPosition_kb.txt"),show_col_types = FALSE)
  
  
  marker_pos_test <- read_table(file.path(data.dir,Dens,"test.map"),
                                col_names = FALSE,show_col_types = FALSE)
  colnames(marker_pos_test) <- c("Chr","locus", "cM", "BPP")
  
  
  # attach the observable position.
  marker_pos_test <- left_join(marker_pos_test, marker_pos[,c("locus","tmp_marker")], by = "locus") %>%
    arrange(BPP)
  
  # To save myself some time we parallelize over ID
  
  # results_list <- for(i in 1:nrow(homo_bin) ) {
  #   result <- process_homo(i, marker_pos_test, homo_bin)
  #   result
  # }
  # 
  # # Combine the results into a single data frame for this generation
  # ID_result <- do.call(rbind, results_list)
  
  # Method 1: Using lapply()
  results_list <- lapply(1:nrow(homo_bin), function(i) {
    process_homo(i, marker_pos_test, homo_bin)
  })
  
  # Combine the results into a single data frame
  ID_result <- do.call(rbind, results_list)
  
  ID_result_2 <- ID_result %>%
    group_by(Generation, locus, tmp_marker ) %>%
    dplyr::summarise(obs_homo = sum(obs_hom),
              n_ind = n() ) %>%
    mutate(Ne = NE,
           Rep = Rep,
           Dens = Dens)
  
  write.table(ID_result_2, file = file.path(data.dir,Dens,"homozyg_test_markers.txt"),sep = "," , quote = FALSE, row.names = FALSE)
}





