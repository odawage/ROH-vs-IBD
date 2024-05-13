
# ====
#   Packages 
# ====
library(data.table)
library(readr)
library(stringi)
library(tidyverse)


# ====

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]

#name <- "NE200_100cm_rep0"

dataset <- name
# Extract if it's a full or half density scenario 
if (grepl("_50", name)) {
  density  <- "Half"
}else{
  density = "Full"
}

# Because of ADAMs file naming convention etc
# I need to pull out replica information 
rep <- as.numeric(gsub("[^0-9]", "", sub(".*_rep", "", dataset)))

#============
# Dir names 
#===============
ibd.dir <- paste0("Dataset/",dataset,"/IBD_segments")

R.dir <- paste0("Dataset/",dataset,"/R_output")

file.dir <- paste0("Dataset/",dataset)


#============
# Make dirs 
#===============

dir.create(ibd.dir)
#dir.create(paste0(ibd.dir,"/Half_Dens"))
#dir.create(paste0(ibd.dir,"/Full_Dens"))

dir.create(R.dir)
#dir.create(paste0(R.dir,"/Half_Dens"))
#dir.create(paste0(R.dir,"/Full_Dens"))



# ====
#  datasets
# ====

# SNP data
marker_QTL <- read_csv(paste0(file.dir,"/MarkerQtlGenotypesRep_rel.res.bz2"), 
                       col_types = cols(Geno = col_character()) )

colnames(marker_QTL) <- c("ID","Genome") 

marker_pos <- read.table(paste0(file.dir,"/markerPositionsRep",rep,".res.bz2"), header = TRUE)

# Founder marker data
marker_templated <- read.delim(paste0(file.dir,"/templateGenotypesRep",rep,".res.bz2"),header = FALSE)

templated_pos <- read.table(paste0(file.dir,"/templateMarkerPositions.res.bz2"), header = TRUE)

# Pedigree
ped <- read_table(paste0(file.dir,"/animalsRep",rep,".res.bz2"), 
                  col_types = cols(birthHerd = col_skip(), 
                                   polyTbv1 = col_skip(), qtlTbv1 = col_skip(), 
                                   residual1 = col_skip()))

# ========
#   this is from the script 1: data_wrangling before loop
# ========

# Data analysis generations 
time <- c(20,40,60,80,100)
rel_ID_short <- ped%>%  
  filter(birth %in% time)


# Rate of inbreeding generations 
time <- c(12:17,18:22,32:37,38:42,52:57,58:62,72:77,78:82,89:100)
rel_ID <- ped%>%  
  filter(birth %in% time)


rel_ID_pos <- which(ped$id %in% rel_ID$id) 


marker_templated <- tibble(marker_templated[rel_ID_pos, ])
colnames(marker_templated) <- "V1"



# =======
#   marker postion cM to bpp
# =======


marker_pos <- marker_pos %>% 
  mutate(position_bpp = position*1000000) %>% 
  mutate(tmp_marker = c(1:nrow(marker_pos)))

templated_pos <- templated_pos %>% 
  mutate(position_bpp = position*1000000) %>% 
  mutate(tmp_marker = c(1:nrow(templated_pos)))


write.table(marker_pos, file= paste0(R.dir,"/markerPosition_kb.txt"), row.names = FALSE, sep = ",", quote = FALSE)  

write.table(templated_pos, file= paste0(R.dir,"/true_IBD_position.txt"), row.names = FALSE, sep = ",", quote = FALSE)  



# ==========
#  QTL markers cleanup   
# =========

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

data <- data.frame(ID = marker_QTL$ID, geno = clean_genos)


write_csv(data, file= paste0(R.dir,"/All_ind_obs_QTL.res.bz2"), quote = c("none") )  

# remove the now irrelevant version 
#file.remove(paste0(file.dir,"/MarkerQtlGenotypesRep_rel.res.bz2"))


# ==========
#  templated markers cleanup   
#     This code is slightly different because of the different formatting for the 
#     Founder markers
#=========



## Remove the weird whitespace at the start of the dataset 
test_2 <- data.frame(clean = str_trim(marker_templated$V1, side = c("left"))) 

## select out the ID column 
test_2$ID <- sapply(strsplit(test_2$clean, " "), `[` , 1)

# rename since it the dissapear when forcing to be dataframe 
colnames(test_2) <- c("geno","ID")
sample_true_IBD <- test_2  # re-save for safety 

# couldn't make apply() work so we brute force it with a loop 
for (i in seq(nrow(test_2))) {
  sample_true_IBD$geno[i] <-  paste0(strsplit(test_2$geno[i], " ")[[1]][-1], collapse = " ")
  #cat(i, "\n")  
}


all_ind_true_ibd <- sample_true_IBD

all_ind_true_ibd <-all_ind_true_ibd %>% 
  relocate(ID) %>% 
  arrange(ID)

all_ind_true_ibd_filtered <- all_ind_true_ibd %>% 
  filter(ID %in% rel_ID_short$id)
## save the dataset as a .txt so I dont have to reload that giant file 
write_csv(all_ind_true_ibd_filtered, file= paste(R.dir,"/All_Ind_true_IBD.res.bz2",sep = ""), quote = c("none"))  
write_csv(all_ind_true_ibd, file= paste(R.dir,"/All_Ind_true_IBD_long.res.bz2",sep = ""), quote = c("none"))  





# empty the environment for posterity and peace of mind

rm(list = setdiff(ls(envir = globalenv()), c("dataset","density" ,"file.dir","ibd.dir","R.dir","res.dir","rel_ID","rel_ID_short")), envir = globalenv())



# ====
# datasets 
# ====

marker_pos <- read.csv(paste0(R.dir,"/markerPosition_kb.txt"))
obs_QTL <- read_csv(paste0(R.dir,"/All_ind_obs_QTL.res.bz2"), 
                           col_types = cols(geno = col_character() ))

templated_IBD <- read.csv(paste0(R.dir,"/All_Ind_true_IBD_long.res.bz2"))
temp_IBD_pos <- read.csv(paste0(R.dir,"/true_IBD_position.txt"))



# ==========
#   Make a binary string of all the places the SNPs are homozygous 
#   
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

write_csv(obs_homo,paste0(R.dir, "/all_ind_observed_homozyg_QTL.res.bz2"), quote = c("none"))     



# ======================
# Find all the founder markers that are IBD 
# and write it out as a binary genome 
# 1 = IBD, 0 = non IBD 
# ===================


ID_vec <- templated_IBD$ID
ID_vec <-sort(ID_vec)

obs_ibd <- data.frame()

for (ID in ID_vec) {
  nmarker <- nrow(temp_IBD_pos)
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
temp <- left_join(result_table_bin, temp_IBD_pos[,3:4], by = c("start" = "tmp_marker")) %>% 
  rename(start_bpp = position_bpp)
info_ADAM_ROH_bin <-  left_join(temp, temp_IBD_pos[,3:4], by = c("stop" = "tmp_marker")) %>% 
  rename(stop_bpp = position_bpp)

info_ADAM_ROH_bin <- info_ADAM_ROH_bin %>%  
  arrange(ID) %>% 
  mutate(NSNP = stop - start)

write.table(info_ADAM_ROH_bin, paste0(ibd.dir, "/info_ADAM_roh_bin.txt"), row.names = FALSE, sep = ",", quote = FALSE)    

file.remove(paste0(R.dir,"/All_Ind_true_IBD_long.res.bz2"))


# =========
#   Translate the IBD segments from founder marker map 
#   to the SNP marker map 
# =========

ID_vec <- templated_IBD %>% 
  filter(ID %in% rel_ID_short$id)

ID_vec <- sort(ID_vec$ID)


all_genome <- data.frame()


for (ID in ID_vec){
  info_ind <- info_ADAM_ROH_bin[info_ADAM_ROH_bin$ID %in% ID, ]
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
        ID = ID,
        genome = genome_str,
        NNA = NNA
      ) 
    } #run   
  }else{
    #cat("no run",ID,"\n")
    } 
  
  all_genome <- rbind(all_genome, genome_id)
  #cat(ID, "\n")
} # ID



write_csv(all_genome,  paste0(ibd.dir, "/extend_bin_ADAM_ROH.res.bz2"), quote = c("none"))       







