
# =============
#
#
#
#
#
# =============
library(stringi)
library(dplyr)
library(readr)
library(data.table)
library(tidyr)

# tester <- function(x){
#   test <- paste0("mjau mjau mjau, vi er ",x*x," katter")
#   return(test)
# }
# 
# tester(c(1:5))
# 
# rm(obj)



# =============
#   loading basic files.  
# =============


RohBinGeno <- function(ROH_inf,marker_pos) {
  # ========
  #   Write a binary genome with ROH state 
  #   1 marker in ROH, 0 as non-ROH 
  # =======
  
  ID_vec <- unique(ROH_inf$ID)
  
  ROH_pos <- data.frame()
  
  for (ID in ID_vec) {
    genome <- as.character( array(data = 0, dim = nrow(marker_pos)))  
    
    ROH_inf_ind <- ROH_inf[ROH_inf$ID %in% ID, ] # select out all ROH for one ID
    ind_run_vec <- c(1:nrow(ROH_inf_ind))

    # Then for all markers that are within a run put a 1 on the genome 
    for (run in ind_run_vec) {
      start <- ROH_inf_ind$POS1[run]
      stop <- ROH_inf_ind$POS2[run]
      #find markers in ROH and write 1s
      rel_marker <- marker_pos$tmp_marker[marker_pos$position_bpp >= start & marker_pos$position_bpp <= stop]
      genome[rel_marker] <- paste(1)
      
    }# RUN
    genome_str <-  paste0(genome, collapse = " ")   
    genome_id <- data.table(
      ID = ID,
      genome = genome_str)
    
    ROH_pos <- rbind(ROH_pos, genome_id)
  } # ID 
  
  return(ROH_pos)
}


# ============
# ROH|IBD  - loop
# =============


RohPower <- function(info_ADAM_ROH_bin, PLINK_ROH_bin, marker_pos, tmp.dir) {
  
  ID_vec <- unique(info_ADAM_ROH_bin$ID)
  ID_plink <- unique(PLINK_ROH_bin$ID)
  
  # In the first generations not all individuals have IBD segments 
  # So here I create "empty" genomes for these IDs to make sure all runs normally 
  leftovers <- ID_vec[!ID_vec %in% ID_plink]
  gen <- (paste0(as.character( array(data = 0, dim = nrow(marker_pos))), collapse = " "))
  len <- length(leftovers)
  
  ext <- data.frame(ID = leftovers, 
                    genome = rep(gen,len))
  
  PLINK_ROH_bin <- rbind(PLINK_ROH_bin, ext)
  
  rm("leftovers","gen","len","ext")
  
  PLINK_ROH_bin <- PLINK_ROH_bin %>% 
    distinct(ID, .keep_all = TRUE)
  
  #starting the loopti-loop
  final <- data.frame( ID = character(0), run = character(0), result = character(0)) 
  
  for (i in ID_vec){
    
    IBD_ind  <- PLINK_ROH_bin %>% 
      filter(ID == i) %>% 
      pull(genome)
    
    info_ind <- info_ADAM_ROH_bin %>% 
      filter(ID == i)
    
    ind_run_vec <- c(1:nrow(info_ind))
    ind_result <- data.frame(matrix(NA, nrow = nrow(info_ind), ncol = 3))
    colnames(ind_result) <- c("ID","ROH","result")
    ind_result$ID <- info_ind$ID
    
    segment <- strsplit(IBD_ind, " ")[[1]] # split genome string into a vector 
    segment <- segment[!grepl("^\\s*$", segment)]  #remove all the empty chunks 
    for(run in ind_run_vec) {
      ind_result$ROH[run] <- run # Store what per ID loop it is 
      
      start <- info_ind$start_bpp[run] 
      stop <- info_ind$stop_bpp[run]
      
      # Since the IBD segments are found with a text string function, 
      # Some of the registered segements have only 1 marker. So first i exclude them
      if( start == stop){cat("no run", i, run ,"\n")
      }else{
        rel_marker <- marker_pos$tmp_marker[marker_pos$position_bpp >= start & marker_pos$position_bpp <= stop]
        
        # The next issue is that for very short runs sometimes we wont have a 
        # SNP marker between start and stop. these areas are then removed 
        if(length(rel_marker) == 0){cat("no markers", i, run ,"\n")
        }else{
          # To be able to compare the results to the raw data. I still write the results 
          # Into a genome string format. While ineficient, it  10/10 for troubleshooting 
          run_res <- vector()
          for (marker in rel_marker) {
            
            if (segment[marker] == 1){
              run_res_t<- paste(1)  #is correct 
            } else {
              run_res_t <- paste(0)
            }# else 
            run_res <- c(run_res, run_res_t)
          } # marker  
          ind_result$result[run] <- paste0(run_res, collapse = " ")
        } # else rel_mark
      } # else start stop
    } # run 
    final <- rbind(final, ind_result)
  } # IDs
  
  # then I examine all the genome strings and cook it down to per. segment numbers
  final$Length_SNP <- (stringi::stri_length(final$result) - stringi::stri_count_fixed(final$result, " "))
  final$sum <- sapply(strsplit(final$result," "), function(x) sum(as.numeric(x)))
  final$pros_right <- final$sum/final$Length_SNP*100
  
  write_csv(final, paste0(tmp.dir, "/ROH_Power.csv.gz"), quote = "none")    
  
}


# ============
# IBD|ROH  - loop
# =============


RohTruePos <- function(ROH_inf, obs_IBD, obs_homo,marker_pos, tmp.dir) {
  
  # There is not a 100% match between what IDs are seen in the two datasets 
  # So here I'm just finding out who's missing and giving them a completely 
  # non-IBD 00000... genome 
  if (nrow(ROH_inf) > 0) {
  leftovers <- obs_homo$ID[!obs_homo$ID %in% obs_IBD$ID]
  gen <- (paste0(as.character( array(data = 0, dim = nrow(marker_pos))), collapse = " "))
  len <- length(leftovers)
  
  ext <- data.frame(ID = leftovers, 
                    genome = rep(gen,len))
  obs_IBD <- rbind(obs_IBD, ext)
  obs_IBD <- obs_IBD %>% 
    distinct(ID, .keep_all = TRUE)
  
  
  run_vec <- 1:nrow(ROH_inf)
  all_run_vec <- nrow(ROH_inf)
  all_runs <- data.frame()
  
  
  # looping it over ROHs
  for (run in run_vec) {
    ID <- ROH_inf$ID[run]
    run_id <- ROH_inf$run[run] # ROH info
    homo_gen <- obs_homo$genome[obs_homo$ID %in% ID] # indi. homozyg. state
    homo_geno <- strsplit(homo_gen, " ")[[1]] 
    
    IBD_gen <- obs_IBD$genome[obs_IBD$ID %in% ID] # Indi. IBD state
    IBD_geno <- strsplit(IBD_gen, " ")[[1]] 
    
    run_inf <- ROH_inf[run, ]
    start <- run_inf$SNP1
    stop <- run_inf$SNP2
    
    rel_markers <-  marker_pos$tmp_marker [marker_pos$locus >= start & marker_pos$locus <= stop]
    
    run_markers <- vector()
    for (marker in rel_markers) {
      # number code:
      #  0 = HET
      #  1 = Non HOMO IBD -> can be the masked markers
      #  2 = IBS / AIS 
      #  3 = IBD 
      #  4 = something's wrong 
      if (homo_geno[marker] == "1" && IBD_geno[marker] == "1"){
        run_mark <- paste0("3")  # it is IBD 
      }else if (homo_geno[marker] == "1" && IBD_geno[marker] == "0"){
        run_mark <- paste0("2")  # it is IBS 
      }else  if (homo_geno[marker] == "0" && IBD_geno[marker] == "1"){
        run_mark <- paste0("1")  # it is Non HOM IBD  
      }else if (homo_geno[marker] == "0" && IBD_geno[marker] == "0"){ 
        run_mark<- paste0("0")  # it is hetero 
      }else{ 
        run_mark<- paste0("4")  # something's WRONG 
      } # else 
      run_markers <- c(run_markers, run_mark)
    }# marker 
    run_results <- data.frame(ID = ID, 
                              ROH = run_id, 
                              results = paste0(run_markers, collapse = " "), 
                              N_SNPs = length(rel_markers))
    all_runs <- rbind(all_runs, run_results)
    #cat(round(run/all_run_vec*100, digits = 2), "% \n")
  } # run 
  
  
  
  # Function to count unique characters in a string
  count_unique <- function(string) {
    counts <- table(strsplit(gsub(" ", "", string), ""))
    counts_df <- as.data.frame(counts)
    colnames(counts_df) <- c("Number", "Count")
    
    # Set all numbers from 0 to 4 with count 0 if not present
    all_numbers <- as.character(0:4)
    counts_df <- counts_df %>%
      complete(Number = all_numbers, fill = list(Count = 0))
    
    return(counts_df)
  }
  
  # Apply the count_unique function to all rows and combine results into a data frame
  combined_result <- all_runs %>%
    group_by(ID, ROH) %>%
    summarise(result = list(count_unique(results))) %>%
    unnest(result) %>%
    pivot_wider(names_from = Number, values_from = Count, values_fill = 0)
  
  # Print the combined result
  # and rename the columns
  all_results <- left_join(all_runs, combined_result, by = c('ID', 'ROH') ) %>% 
    rename(Hetero = "0",
           N_HOM_IBD = "1", 
           IBS = "2", 
           IBD = "3",
           X4 = "4")
  
  write_csv(all_results,paste(tmp.dir, "/PLINK_ROH_IBD_vs_IBS.res.bz2", sep=""), quote = c("none"))
  
  }else {
    all_results <- data.frame(
      ID = NA,
      ROH = NA,
      results = NA,
      N_SNPs = NA,
      Hetero = NA,
      N_HOM_IBD = NA,
      IBS = NA,
      IBD = NA,
      X4 = NA  # Changed '4' to 'X4' to avoid issues
    )
    write_csv(all_results,paste(tmp.dir, "/PLINK_ROH_IBD_vs_IBS.res.bz2", sep=""), quote = c("none"))
    
    } # else
  }






