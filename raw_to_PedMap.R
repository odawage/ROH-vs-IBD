
# ===============
#   Making PLINK runable from Rstudio 
# ===============

source("/usr/share/lmod/lmod/init/R")
Sys.setenv(MODULEPATH = '/cluster/modules/all')
module("load PLINK/1.9b_6.17-x86_64")

# =======================
#       Libs 
# ========================
library(tidyverse)
library(readr)

# =========
#   directories 
# =========

NEs = c("NE100","NE200")
Reps = c(1:10)

proj.dir <- paste0("/mnt/project/SimData/Paper_1")


for (Ne in NEs) {
  for (Rep in Reps) {
    dir.create(paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens"))
    dir.create(paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens"))
    # =========
    #   Pedigree and relevant individuals
    # =========
    # load the pedigree
    rep_alt = 1 # reps run i different arrays. files all named 1 
    animalsRep1_res <- read_table( paste0(proj.dir,"/",Ne,"/Rep",Rep,"/animalsRep",rep_alt,".res.bz2"),
                                   col_types = cols(id = col_character(),
                                                    sire = col_character(), dam = col_character(),
                                                    X11 = col_skip()))
    # Time 1: animals used for most of the plots and tables filter relevant individuals
    rel_ID_1 <- animalsRep1_res %>%
      filter(birth %in% c(10,20,50,100))
    # Time 2: animals used for rate of inbreedng
    rel_ID_2 <- animalsRep1_res %>%
      filter(birth %in% c(0:30,40:60,80:100))
    # write out the ped of only the Time 1 individuals
    write.table(rel_ID_1, file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/pedigree_rel_generations.txt"), row.names = FALSE, quote = FALSE, sep = " " )
   
     # =========
    #   SNP data
    # =========
    MarkerQtl <- read_table(paste0(proj.dir,"/",Ne,"/Rep",Rep,"/MarkerQtlGenotypesRep",rep_alt,".res.bz2"),
                            col_names = FALSE,
                            col_types = cols(X1 = col_character(),
                                             X2 = col_character()))
    colnames(MarkerQtl) <- c("ID","Geno") # Rename
    #filter to save some time and space
    MarkerQtl_rel <- MarkerQtl %>%
      filter(ID %in% rel_ID_1$id)
    MarkerQtl_rel_2 <- MarkerQtl %>%
      filter(ID %in% rel_ID_2$id)
    write_csv(MarkerQtl_rel,paste0(proj.dir,"/",Ne,"/Rep",Rep,"/MarkerQtlGenotypesRep_rel.res.bz2"), quote = c("none") )
    rm(MarkerQtl)
    
    # =========
    #   SNP positions
    # =========
    markerPositionsRep1 <- read_table(paste0(proj.dir,"/",Ne,"/Rep",Rep,"/markerPositionsRep",rep_alt,".res.bz2",sep = "") ) %>%
      relocate(chromosome)
    data_pos <- markerPositionsRep1 %>%
      mutate(BPP = position*1000000)
    rm(markerPositionsRep1)
   
     # ===========
    #   Full map, 50_Full_map
    # ===========
    set.seed(11)
    data_pos_50 <- data_pos %>%
      sample_frac(0.5) %>%
      arrange(locus)
    ID_vec_50 <- sort(as.numeric(data_pos_50$locus))
    write_csv(data_pos, file =  paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens/markerPositionsRep1.res.bz2"), quote = c("none")  )
    write.table(data_pos, file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens/full.map"),
                row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE )
    
    write.table(data_pos_50, file =  paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/full.map"),
                row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE )
    write_csv(data_pos_50,file =  paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/markerPositionsRep1.res.bz2"), quote = c("none")  )
   
     # =============
    #   Prepare for MAF filtering for gen 100
    # ============
    #select individuals in the last generation
    ID_max_gen <- animalsRep1_res %>%
      filter(birth == max(animalsRep1_res$birth)) %>%
      select(id)
    # SNP id and observable position is not the same. so i make a relative position column
    # data_pos
    data_pos <- data_pos %>%
      mutate(tmp_pos = c(1:nrow(data_pos)))
    # data_pos_50
    data_pos_50 <- data_pos_50 %>%
      mutate(tmp_pos = c(1:nrow(data_pos_50)))
   
     # =========
    # Here I make the first information chunk needed in the .ped file
    info <- data.frame(matrix(nrow = nrow(MarkerQtl_rel), ncol = 6, 0))
    info[, 2] <- as.numeric(MarkerQtl_rel$ID)
    info_max <- info %>%
      filter(X2 %in% ID_max_gen$id)
    info_rel_gen <- info %>%
      filter(X2 %in% rel_ID_1$id)
    #find where the coordiantes for where the individuals are placed
    position_max <- which(info[ ,2] %in% ID_max_gen$id)
    position_rel_gen <- which(info[ ,2] %in% rel_ID_1$id)
    # info time 1
    info <- unite(info, info, c(1:6), sep = " ")
    # info max gen
    info_max <- unite(info_max, info_max, c(1:6), sep = " ")
    all_SNP <- c(1:(str_length(MarkerQtl_rel$Geno[1])/2))
    #Find all markers without positional data
    leftovers <- all_SNP[-data_pos$locus]
    # positions without positionaldata as haplotype position
    no_pos <- sort( c(leftovers*2-1, leftovers*2))
    # MarkerQTL for max gen
    data <- MarkerQtl_rel %>%
      filter(ID %in% ID_max_gen$id)
   
     # loop over all individuals and remove the markers that don't have positional data
    for (i in seq_len(nrow(data))) {
      seq <- data$Geno[i]
      data[i,2] <- paste(unlist(strsplit(seq,""))[-no_pos], collapse = "")
    }
    
    gen <- data
    rm(data)
    gen$Geno <- gen$Geno%>%  # add 0 to the head
      str_replace_all("(.{1})", "\\1 ") %>% # insert spaces after every to digits
      str_trim() # remove the space in the end if it's there
    info_max$Geno <- gen$Geno
    write.table(info_max, file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/full_gen100.ped"),row.names = FALSE, col.names = FALSE,sep = " ",quote = FALSE )
   
     rm("gen", "info_max")
    
    # MarkerQTL for all rel gen
    data_rel <- MarkerQtl_rel
    for (i in seq_len(nrow(data_rel))) {
      seq <- data_rel$Geno[i]
      data_rel[i,2] <- paste(unlist(strsplit(seq,""))[-no_pos], collapse = "")
    }
    
    rm("seq")
    
    gen_rel <- data_rel
    gen_rel$Geno <- gen_rel$Geno%>%  # add 0 to the head
      str_replace_all("(.{1})", "\\1 ") %>% # insert spaces after every to digits
      str_trim() # remove the space in the end if it's there
    info_rel_gen$Geno <- gen_rel$Geno
    write.table(info_rel_gen, file=paste0(proj.dir,"/",Ne,"/Rep",Rep,"/full_rel_gen.ped"),row.names = FALSE, col.names = FALSE,sep = " ",quote = FALSE )
    
    rm("gen_rel","data_rel", "info_rel_gen")
    
    # ===============
    #   PLINK - MAF
    # ==============
    
    plink.dir <-paste0(proj.dir,"/",Ne,"/Rep",Rep)
    system(paste0("plink --map ",plink.dir,"/Full_Dens/full.map --ped ",plink.dir,"/full_gen100.ped --maf 0.05 --make-bed --recode --out ",plink.dir,"/gen_100_maf"))
    
    # =========
    #   Then we load in the dataset that tells us what SNPS that can be used for the subset
    # ===========
    
    map_maf <- read_table(paste0(plink.dir,"/gen_100_maf.map"),
                          col_names = FALSE)
    colnames(map_maf) <- c("chr","locus", "CM","BPP" )
    # delete the output MAF files
    
    file.remove(dir(
      plink.dir,
      pattern = "^gen_100_",
      full.names = TRUE
    ))
    
    # =========
    #   Then we select 30% of the MAFed SNPs to be test markers
    #   And write ot the test and train .map files
    # =========
    
    # sample of test marker for the full dens
    set.seed(151)
    test <- map_maf %>%
      sample_frac(0.3) %>%
      arrange(locus)
    # remove the test snps
    test_ID <- test$locus
    train <- data_pos[!data_pos$locus %in% test_ID, ]
    # write out
    write.table(train[,1:4], file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens/train.map"), row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE )
    write.table(train, file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens/train_temp.map"), row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE )
    write.table(test, file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens/test.map"), row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE )
    
    # sample out the test markers for the low density scenario
    # ========================================================
    
    set.seed(151)
    test_50 <- map_maf%>%
      filter(locus %in% ID_vec_50) %>%
      sample_frac(0.3) %>%
      arrange(locus)
   
     # Remove the test markers
    test_ID_50 <- test_50$locus
    train_50 <- data_pos_50[!data_pos_50$locus %in% test_ID_50, ]
    write.table(train_50[,1:4], file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/train.map"), row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE )
    write.table(train_50, file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/train_tmp.map", sep = ""), row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE )
    write.table(test_50, file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/test.map", sep = ""), row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE )
   
     # =============
    #   Then we make three .ped files
    # =============
   
     #  one with all marker
    #  one with the train markers
    #  one with the test markers
    # again column with observable position
    data_pos <- data_pos %>%
      mutate(tmp_pos = c(1:nrow(data_pos)))
    data_pos_50 <- data_pos_50 %>%
      mutate(tmp_pos = c(1:nrow(data_pos_50)))
    str_length(MarkerQtl_rel_2$Geno[1])
    # find what markers that are missing positional data,
    # And write their observable position into no_pos
    all_SNP <- c(1:(str_length(MarkerQtl_rel_2$Geno[1])/2))
    leftovers <- all_SNP[-data_pos$locus]
    no_pos <- sort( c(leftovers*2-1, leftovers*2))
    # To find the rate of inbreeding i need more generations
    # So here i use the Time_2 generations to select for the .ped
    data <- MarkerQtl_rel_2
    data <- data %>%
      filter(ID %in% rel_ID_2$id)
    for (i in seq_len(nrow(data))) {
      seq <- data$Geno[i]
      data[i,2] <- paste(unlist(strsplit(seq,""))[-no_pos], collapse = "")
    }
    data_full <- data
    nchar(MarkerQtl_rel_2$Geno[1])
    nchar(data$Geno[1])
   
     # ================
    # Then we need to subset that shit.
    # Here we use the new marker pos label (tmp_pos) since it looks at the observable position
    # in the new cleaned up genome.
    # the test markers are saved as its own text string, the [dataset]_test.map file
    # can then be used to translate between test relative position to the genome position.
    
    data_2 <- data
    
    # =====
    #   50k ped subset
    # =====
    
    K50_pos <- data_pos$tmp_pos[data_pos$locus %in% ID_vec_50]
    K50_pos <- sort( c(K50_pos*2-1, K50_pos*2))
    data_50k <- data.frame(ID = data$ID, Geno = 0 )
    
    for (i in seq_len(nrow(data))) {
      seq <- data$Geno[i]
      data_50k[i, 1] <- data$ID[i]
      data_50k[i,2] <- paste(unlist(strsplit(seq,""))[K50_pos], collapse = "")
    }
    
    nchar(data$Geno[1])
    nchar(data_50k$Geno[1])
    write_csv(data_50k, file = paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/MarkerQtlGenotypesRep1.res.bz2") , quote = c("none")  )
    
    # ========
    # then subset test for the full data
    # ========
    
    test_pos <- data_pos$tmp_pos[data_pos$locus %in% test_ID]
    test_pos <- sort( c(test_pos*2-1, test_pos*2))
    data_test <- data.frame(ID = data$ID, Geno = 0 )
    for (i in seq_len(nrow(data))) {
      seq <- data$Geno[i]
      data_test[i, 2] <- paste(unlist(strsplit(seq,""))[test_pos], collapse = "")
      data_test[i, 1] <- i
      data[i,2] <- paste(unlist(strsplit(seq,""))[-test_pos], collapse = "")
      print(nchar(data$Geno[i]))
    }
    
    nchar(data$Geno[1])
    nchar(data_test$Geno[1])
    nchar(data$Geno[1]) + nchar(data_test$Geno[1])
    
    # ========
    # then subset test for the 50K data
    # ========
    
    test_pos_50k <- data_pos_50$tmp_pos[data_pos_50$locus %in% test_ID_50]
    test_pos_50k <- sort( c(test_pos_50k*2-1, test_pos_50k*2))
    data_test_50k <- data.frame(ID = data_50k$ID, Geno = 0 )
    for (i in seq_len(nrow(data))) {
      seq <- data_50k$Geno[i]
      data_test_50k[i, 2] <- paste(unlist(strsplit(seq,""))[test_pos_50k], collapse = "")
      data_test_50k[i, 1] <- i
      data_50k[i,2] <- paste(unlist(strsplit(seq,""))[-test_pos_50k], collapse = "")
    }
    nchar(data_50k$Geno[1])
    nchar(data_test_50k$Geno[1])
    nchar(data_50k$Geno[1]) + nchar(data_test_50k$Geno[1])
    
    # ===============
    # To make a readable as a .ped the format needs to be
    # 12 22 11 21, so we need to add the spaces
    # ===============
    
    gen_full <- data_full
    gen_full$Geno <- gen_full$Geno%>%  # add 0 to the head
      str_replace_all("(.{1})", "\\1 ") %>% # insert spaces after every two digits
      str_trim() # remove the space in the end if it's there
    gen <- data
    gen$Geno <- gen$Geno%>%  # add 0 to the head
      str_replace_all("(.{1})", "\\1 ") %>% # insert spaces after every two digits
      str_trim() # remove the space in the end if it's there
    gen_50k <- data_50k
    gen_50k$Geno <- gen_50k$Geno%>%  # add 0 to the head
      str_replace_all("(.{1})", "\\1 ") %>% # insert spaces after every two digits
      str_trim() # remove the space in the end if it's there
    gen_test <- data_test
    gen_test$Geno <- gen_test$Geno%>%  # add 0 to the head
      str_replace_all("(.{1})", "\\1 ") %>% # insert spaces after every two digits
      str_trim() # remove the space in the end if it's there
    gen_test_50k <- data_test_50k
    gen_test_50k$Geno <- gen_test_50k$Geno%>%  # add 0 to the head
      str_replace_all("(.{1})", "\\1 ") %>% # insert spaces after every two digits
      str_trim() # remove the space in the end if it's there
    
    # Then the info columns at the start
    # =================
    info <- data.frame(matrix(nrow = nrow(data), ncol = 6, 0))
    info[ , 2] <- as.numeric(data$ID)
    info_split <- tibble(info$X2)
    info <- unite(info, info, c(1:6), sep = " ")
    info_split$Geno <- gen_full$Geno
    # use the same id info for all the output files
    info_full <- info
    info_test <- info
    info_50k <- info
    info_test_50k <- info
    # stick in the genome data
    info_full$Geno <- gen_full$Geno
    info$Geno <- gen$Geno
    info_test$Geno <- gen_test$Geno
    info_50k$Geno <- gen_50k$Geno
    info_test_50k$Geno <- gen_test_50k$Geno
    # write it all out
    write.table(info_full, file= paste0(proj.dir,"/",Ne,"/Rep",Rep,"/full.ped"), row.names = FALSE, col.names = FALSE,sep = " ",quote = FALSE )
    write.table(info_split, file =  paste0(proj.dir,"/",Ne,"/Rep",Rep,"/full_split.ped"), row.names = FALSE, col.names = FALSE,sep = ",",quote = FALSE)
    write.table(info, file =  paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens/train.ped"), row.names = FALSE, col.names = FALSE,sep = " ",quote = FALSE )
    write.table(info_test, file =  paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens/test.ped"), row.names = FALSE, col.names = FALSE,sep = " ",quote = FALSE)
    write.table(info_50k, file =  paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/train.ped"),row.names = FALSE, col.names = FALSE,sep = " ",quote = FALSE )
    write.table(info_test_50k, file =  paste0(proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/test.ped"), row.names = FALSE, col.names = FALSE,sep = " ",quote = FALSE)
   
     # ======================
    #   Then PLINK again to make  the beforeQC .bim/.bed/.fam files
    #   that's getting sent to the --homozyg function
    # ======================
   
     system(paste0("plink --file " ,proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens/train --make-bed --out " ,proj.dir,"/",Ne,"/Rep",Rep,"/Full_Dens/beforeQC"))
    system(paste0("plink --file ",proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/train --make-bed --out " ,proj.dir,"/",Ne,"/Rep",Rep,"/Half_Dens/beforeQC"))
  } # rep
} # NE

