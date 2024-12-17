# ===========================
# Calling PLINK from Rstudio

# source("/usr/share/lmod/lmod/init/R")
# Sys.setenv(MODULEPATH = '/cluster/modules/all')
# module("load PLINK/1.9b_6.17-x86_64")
# plink = "plink"

# =======================
# Libraries 


library(data.table)
library(plyr)
library(openxlsx)
library(ggplot2)
library(readr)
# 
# input.dir <- "/mnt/users/odwa/PLINK/Code/Dataset"
# dataset = "NE100_100cm_rep1"
# inped = "beforeQC_50"
# tmp.dir = "/net/fs-2/scale/OrionStore/Scratch/odwa/paper-1"

# Gorssen(input.dir, dataset,inped,tmp.dir )


Gorssen <- function(input.dir, dataset, inped, tmp.dir, dens, err ) {
  plink = "plink"
  # =
  # 1. Window settings:
  #   - window length, scores' SNPS to decide if in segment or not
  # =

  # fixed para
  het <-0  #--homozyg-window-het AND --homozyg-het 
  mis <-0  #--homozyg-window-missing
  Nout <- 2 #The number of outer SNPs you don't want to allow for ROH-analysis (calculation of windowthreshold)
  # Adjusted
  windowsnp <- 0 ##--homozyg-window-snp => this parameter will be adjusted within script to parameter l!
  windowthreshold <- 0.05 #--homozyg-window-threshold => this parameter will be adjusted within script (t=floor(Nout+1L,3)!

  #=
  #2.Segment settings: - #After segment calling via window settings: extra requirements for segments!
  #=

  # fixed
  kb <- 1000 #--homozyg-kb: the minimum segment length in kb
  gap <- 500 #--homozyg-gap: the maximal gap in kb between two SNPs within an ROH
  density <- 60  #--homozyg-density #the minimal average density of an ROH (expressed as 1SNP/xxx kb, so in our example: minimal 1SNP/150kb)
  # Adjusted
  l<-0    #--homozyg-snp => this parameter will be adjusted within script!

  #====================================
  # dataset
  #====================================

  name <- dataset
  #name = "NE100_100cm_rep1"
  if (grepl("_50", name)) {
    data_name  <- sub("_50", "", name)
  }else{
    data_name = name
  }

  #=====
  #  Create directorie 
  # =======

  dir.create(paste0(tmp.dir,"/",data_name))

  if (grepl("_50", dens)) {
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens/Meyermans"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens/Meyermans/",err))
    
    plink.out <- paste0(tmp.dir,"/",data_name,"/Half_Dens/Meyermans/",err)  

  }else{
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens/Meyermans/"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens/Meyermans/",err))

    plink.out <- paste0(tmp.dir,"/",data_name,"/Full_Dens/Meyermans/",err)
  }


  # ========================
  #   The parts needed for the formulas 
  # ========================

  fam <- read.table(paste0(input.dir,"/",inped,".fam"),h=F)

  #load bim file: determine number of chromosomes

  map <-read.table(paste0(input.dir,"/",inped,".bim"),h=F)
  


  system(paste0(plink," --allow-extra-chr --bfile ",input.dir,"/",inped," --hardy --make-bed --out ",plink.out,"/ROH_analyse"))

  #Load results
  plink.hwe <- read.table(file=paste0(plink.out,"/ROH_analyse.hwe"),fill=TRUE,header=TRUE)
  #write  observed heterozygosity
  heter <- mean(plink.hwe$O.HET.,na.rm = T) #

  #===
  # calculate 
  #===

  #l=log(0.05/(number_of_SNPs*number of animals))/log(1-mean_heterozygosity)
  l <- (log(0.05/(nrow(map)*nrow(fam)))/log(1-heter))

  #Round l parameter to an integer number
  l <- round(l)
  if (is.na(l)){l <- 50} #In a rare occasion l is set to NA, set l = 50 in this case  
  #Change window-snp to l
  windowsnp <- l ##--homozyg-window-snp

  #Calculate windowthreshold based on formula from Meyermans&Gorssen et al, 2020:
  windowthreshold<-floor(1000*((Nout+1)/l))/1000



  #===
  #4. ROH analysis in Plink (--homozyg)
  #===

  command <- paste0(plink," --allow-extra-chr --bfile ",plink.out,"/ROH_analyse", "  --homozyg-window-het ",het,
  " --homozyg-het ",het," --homozyg-window-missing ",mis," --homozyg-snp ",l," --homozyg-kb ",kb," --homozyg-window-snp ",
  windowsnp,"  --homozyg-window-threshold ",windowthreshold,"  --homozyg-density ",density," --homozyg-gap ",gap," --out ",plink.out,"/ROH_analyse")
  
  readr::write_lines(command, file = paste0(plink.out,"/PLINK_mey_commands.txt") )

  system(command)
}


Norm <- function(input.dir, dataset,inped,tmp.dir, dens,err) {
  #====================================
  # dataset
  #====================================
  
  data_name <- dataset
  #name = "NE100_100cm_rep1"
  # if (grepl("_50", name)) {
  #   data_name  <- sub("_50", "", name)
  # }else{
  #   data_name = name
  # }
  
  #=====
  #  Create directorie 
  # =======
  
  dir.create(paste0(tmp.dir,"/",data_name))
  
  if (grepl("_50", dens)) {
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens/Norm"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens/Norm/",err))
    
    plink.out <- paste0(tmp.dir,"/",data_name,"/Half_Dens/Norm/",err)  
    
  }else{
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens/Norm"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens/Norm/",err))
    
    plink.out <- paste0(tmp.dir,"/",data_name,"/Full_Dens/Norm/",err)   
  }
  system(paste0("plink --bfile ",input.dir,"/",inped, " --homozyg --homozyg-snp 50 --homozyg-kb 500 --homozyg-density 100 --homozyg-gap 1000 --homozyg-het 0  --homozyg-window-snp 50 --homozyg-window-het 0 --homozyg-window-missing 1 --homozyg-window-threshold 0.05 --out ",plink.out, "/ROH_analyse"))
  
}



Norm_small <- function(input.dir, dataset,inped,tmp.dir, dens,err) {
  data_name <- dataset

  
  #=====
  #  Create directorie 
  # =======
  
  dir.create(paste0(tmp.dir,"/",data_name))
  
  if (grepl("_50", dens)) {
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens/Norm_small"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens/Norm_small/",err))
    
    plink.out <- paste0(tmp.dir,"/",data_name,"/Half_Dens/Norm_small/",err)  
    
  }else{
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens/Norm_small"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens/Norm_small/",err))
    
    plink.out <- paste0(tmp.dir,"/",data_name,"/Full_Dens/Norm_small/",err)   
  }
  system(paste0("plink --bfile ", input.dir,"/",inped," --homozyg --homozyg-snp 20 --homozyg-kb 500 --homozyg-density 100 --homozyg-gap 1000 --homozyg-het 0 --homozyg-window-snp 20 --homozyg-window-het 0 --homozyg-window-missing 1 --homozyg-window-threshold 0.05 --out ", plink.out,"/ROH_analyse"))
  
  }



Default <- function(input.dir, dataset, inped,tmp.dir, dens,err) {
  name <- dataset
  #name = "NE100_100cm_rep1"
  if (grepl("_50", dens)) {
    data_name  <- sub("_50", "", name)
    
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens/Default"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Half_Dens/Default/",err))
    
    plink.out <- paste0(tmp.dir,"/",data_name,"/Half_Dens/Default/",err)
  }else{
    data_name = name
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens/Default"))
    dir.create(paste0(tmp.dir,"/",data_name,"/Full_Dens/Default/",err))
    
    
    plink.out <- paste0(tmp.dir,"/",data_name,"/Full_Dens/Default/",err) 

  }
  system(paste0("plink --bfile ", input.dir,"/",inped ," --homozyg --out ",plink.out ,"/ROH_analyse"))
         
}
  

