
# ====
#   Packages 
# ====
library(data.table)
library(readr)
library(stringi)
library(tidyverse)

# ====


# ====
#  Pathing 
# ===


args <- commandArgs(trailingOnly = TRUE)
path <- args[1]
NE <- args[2]
Rep <- args[3]
# 
# path <- "/mnt/project/SimData/Paper_1"
# NE <- "NE100"
# Rep <- 1

data.dir <- paste0(path,"/",NE,"/Rep",Rep)



# =======
#   marker postion cM to bpp
# =======

marker_pos <- read.table(paste0(data.dir,"/markerPositionsRep1.res.bz2"), header = TRUE)
templated_pos <- read.table(paste0(data.dir,"/templateMarkerPositions.res.bz2"), header = TRUE)


marker_pos <- marker_pos %>% 
  mutate(position_bpp = position*1000000) %>% 
  mutate(tmp_marker = c(1:nrow(marker_pos)))

templated_pos <- templated_pos %>% 
  mutate(position_bpp = position*1000000) %>% 
  mutate(tmp_marker = c(1:nrow(templated_pos)))


write.table(marker_pos, file= paste0(data.dir,"/markerPosition_kb.txt"), row.names = FALSE, sep = ",", quote = FALSE)

write.table(templated_pos, file= paste0(data.dir,"/true_IBD_position.txt"), row.names = FALSE, sep = ",", quote = FALSE)


# ==========
#  QTL markers cleanup
# =========
# SNP data
marker_QTL <- read_csv(paste0(data.dir,"/MarkerQtlGenotypesRep_rel.res.bz2"),
                       col_types = cols(Geno = col_character()) )
colnames(marker_QTL) <- c("ID","Genome")
# I want to use the MarkerQtlGenotypesRep file format,
# but I still need to remove the markers without positions

# So find the leftovers
n_QTLs <- nchar(marker_QTL$Genome[1])/2
leftovers <- c(1:n_QTLs)[-marker_pos$locus]

clean_genos <- vector()

for (i in 1:nrow(marker_QTL)) {
  # split the genome every second marker
  split_geno <- strsplit(marker_QTL$Genome[i], "(?<=.{2})", perl = TRUE)[[1]]
  if (length(leftovers)!= 0) {

    gened <- split_geno[-leftovers]
  } else{
    gened <- split_geno
  }
  clean_geno <- paste0(gened, collapse = "")
  clean_genos <- c(clean_genos, clean_geno)
  #cat(i, "\n")
}

obs_QTL <- data.frame(ID = marker_QTL$ID, geno = clean_genos)
write_csv(obs_QTL, file= paste0(data.dir,"/All_ind_obs_QTL.res.bz2"), quote = c("none") )

cat("Done with QTL cleanup")
rm(marker_QTL)

# ==========
#  templated markers cleanup
#     This code is slightly different because of the different formatting for the
#     Founder markers
# =========
marker_templated <- read.delim(paste0(data.dir,"/templateGenotypesRep1.res.bz2"),header = FALSE)
rel_ped <- fread(paste0(data.dir,"/pedigree_rel_generations.txt"))

## Remove the weird whitespace at the start of the dataset
sample_true_IBD <- data.frame(clean = str_trim(marker_templated$V1, side = c("left")))

## select out the ID column
sample_true_IBD$ID <- sapply(strsplit(sample_true_IBD$clean, " "), `[` , 1)

# rename since it the dissapear when forcing to be dataframe
colnames(sample_true_IBD) <- c("geno","ID")

# couldn't make apply() work so we brute force it with a loop
for (i in seq(nrow(sample_true_IBD))) {
  sample_true_IBD$geno[i] <-  paste0(strsplit(sample_true_IBD$geno[i], " ")[[1]][-1], collapse = " ")
  #cat(i, "\n")
}

sample_true_IBD <- sample_true_IBD %>%
  relocate(ID) %>%
  arrange(ID)

all_ind_true_ibd_filtered <- sample_true_IBD %>%
  filter(ID %in% rel_ped$id)
## save the dataset as a .txt so I dont have to reload that giant file
write_csv(all_ind_true_ibd_filtered, file= paste(data.dir,"/All_Ind_true_IBD.res.bz2",sep = ""), quote = c("none"))

#write_csv(sample_true_IBD, file= paste(data.dir,"/All_Ind_true_IBD_long.res.bz2",sep = ""), quote = c("none"))

templated_IBD <- sample_true_IBD

rm(list = setdiff(ls(envir = globalenv()), c("NE","Rep","data.dir","rel_ped","templated_IBD","marker_pos","templated_pos","obs_QTL")), envir = globalenv())

cat("Done with founder marker cleanup")


# ==========
#   Make a binary string of all the places the SNPs are homozygous
# ==========

ID_vec <- obs_QTL$ID
obs_homo <- data.frame()

for (ID in ID_vec) {
  # Empty dataset and vector
  genome <-  as.character( array(data = 0, dim = nrow(marker_pos)))
  NSNP <- 0

  # Extract genome and split into vector
  geno_temp <- obs_QTL[obs_QTL$ID %in% ID, ]
  geno <- strsplit(geno_temp$geno, "")[[1]]

  for (marker in 1:nrow(marker_pos)) {
    # If the marker is homozyg. write 1, else write 0
    if (geno[marker*2-1] == geno[marker*2]){
      genome[marker] <- 1
      NSNP <- NSNP +1
    } else {
      genome[marker] <- 0
    } # else
  } # marker
  NNA <- sum(is.na(genome))
  genome_str <- paste0(genome, collapse = " ")
  genome_ind <- data.frame(
    ID = ID,
    genome = genome_str,
    NSNP = NSNP,
    NNA = NNA
  )
  obs_homo <- rbind(obs_homo, genome_ind)
  #cat("we're on ind", ID, sep = "\n")

}# ID

write_csv(obs_homo,paste0(data.dir,"/all_ind_observed_homozyg_QTL.res.bz2"), quote = c("none"))

cat("Done with binary homozyg")

# ======================
# Find all the founder markers that are IBD
# and write it out as a binary genome
# 1 = IBD, 0 = non IBD
# ===================


ID_vec <- templated_IBD$ID
ID_vec <-sort(ID_vec)

obs_ibd <- data.frame()

for (ID in ID_vec) {
  nmarker <- nrow(templated_pos)
  genome <-  as.character( array(data = 0, dim = nmarker))

  geno_temp <- templated_IBD[templated_IBD$ID %in% ID, ]
  geno <- strsplit(geno_temp$geno, " ")[[1]] #split
  geno <- geno[!grepl("^\\s*$", geno)] # remove empty elements
  #genome <- character(nmarker)
  genome_bin <- character(nmarker)
  NSNP <- 0

  for (marker in 1:nmarker) {
    if (geno[marker*2 - 1] == geno[marker*2]) {
      #genome[marker] <- geno[marker*2]
      genome_bin[marker] <- 1
      NSNP <- NSNP +1
    }else{
      #genome[marker] <- 0
      genome_bin[marker] <- 0
    }#else
  }# marker
  NNA <- sum(is.na(genome)) # measures how often something is wrong
  #genome_str <- paste0(genome, collapse = " ")
  genome_bin_str <- paste0(genome_bin, collapse = " ")

  genome_ind <- data.frame(
    ID = ID,
    #genome = genome_str,
    genome_bin = genome_bin_str,
    NSNP = NSNP,
    NNA = NNA
  )
  obs_ibd <- rbind(obs_ibd, genome_ind)
  # cat("we're on ind", ID, sep = "\n")
}

# give some individual stats
obs_ibd$Length_SNP <- (stringi::stri_length(obs_ibd$genome_bin) - stringi::stri_count_fixed(obs_ibd$genome_bin, " "))
obs_ibd$sum <- sapply(strsplit(obs_ibd$genome_bin," "), function(x) sum(as.numeric(x)))
obs_ibd$pros_IBD <- obs_ibd$sum/obs_ibd$Length_SNP*100

cat("Done with binary founder")


# ===============
# Find the start and stop of the binary IBD segments
# ===========

# Create an empty data.table to store the results
result_table_bin <- data.frame()

# Loop over each Individual
for (i in seq_along(obs_ibd$genome_bin)) {
  geno <- obs_ibd$genome_bin[i]
  geno <- str_remove_all(geno, "\\s+")
  # Find all positions of the continuous "1" areas in the genome
  matches <- str_locate_all(geno, "1+")
  matches <- matches[[1]]

  # Check if there are any areas of continuous "1"
  if (nrow(matches) > 0) {
    # Create a data.table with the results for the current individual
    array_result <- data.frame(
      ID = obs_ibd$ID[i],
      run = 1:nrow(matches),
      start = matches[, "start"],
      stop = matches[, "end"]
    )

    # Append the results for the current array to the result table
    result_table_bin <- rbind(result_table_bin, array_result)
  } else {
    cat("No ROH found in individual", i, "\n")
  }
}

# since start and stop is only based on observable position we have to
# attach the positional data using the tmp_markers from earlier
temp <- left_join(result_table_bin, templated_pos[,3:4], by = c("start" = "tmp_marker")) %>%
  rename(start_bpp = position_bpp)
info_ADAM_ROH_bin <-  left_join(temp, templated_pos[,3:4], by = c("stop" = "tmp_marker")) %>%
  rename(stop_bpp = position_bpp)

info_ADAM_ROH_bin <- info_ADAM_ROH_bin %>%
  arrange(ID) %>%
  mutate(NSNP = stop - start)

#write.table(info_ADAM_ROH_bin, paste0(data.dir,"/info_ADAM_roh_bin.txt"), row.names = FALSE, sep = ",", quote = FALSE)
write.table(info_ADAM_ROH_bin, file.path(data.dir,"info_ADAM_roh_bin.txt"), row.names = FALSE, sep = ",", quote = FALSE)

#file.remove(paste0(data.dir,"/All_Ind_true_IBD_long.res.bz2"))

cat("done with ADAM segment coords")




# =========
#   Translate the IBD segments from founder marker map 
#   to the SNP marker map 
# =========

info_ADAM_ROH_bin <- fread( file.path(data.dir,"info_ADAM_roh_bin.txt"))

ID_vec <- rel_ped$id

all_genome <- data.frame()


for (i in ID_vec){
  info_ind <- filter(info_ADAM_ROH_bin, ID == i)
  ind_run_vec <- c(1:nrow(info_ind))

  # make the empty genome 
  genome <- as.character( array(data = 0, dim = nrow(marker_pos)))
  if (length( ind_run_vec) > 0){
    for (run in ind_run_vec) {
      start <- info_ind$start_bpp[run]
      stop <-  info_ind$stop_bpp[run]
      # find all the markers between start and stop of a IBD segment
      rel_marker <- marker_pos$tmp_marker [marker_pos$position_bpp >= start & marker_pos$position_bpp <= stop]
      for (marker in rel_marker){
        genome[marker] <- paste(1)
        #genome[is.na(genome)] <- 0
      } #marker 
      NNA <- sum(is.na(genome)) #check if some spots have problems 
      genome_str <-  paste0(genome, collapse = " ")   
      
      genome_id <- data.table(
        ID = i,
        genome = genome_str,
        NNA = NNA
      ) 
    } #run   
  }else{
    cat("no run",i,"\n")
  } 
  
  all_genome <- rbind(all_genome, genome_id)
  #cat(i, "\n")
} # ID



write_csv(all_genome,  file.path(data.dir, "extend_bin_ADAM_ROH.res.bz2"), quote = c("none"))       

