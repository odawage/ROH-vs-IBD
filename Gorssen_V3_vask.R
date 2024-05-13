
# 
# source("/usr/share/lmod/lmod/init/R")
# Sys.setenv(MODULEPATH = '/cluster/modules/all')
# module("load PLINK/1.9b_6.17-x86_64")



plink = "plink"
#===
#Load and install following functions and packages
#===

#install.packages('data.table')
library(data.table)
#install.packages('plyr')
library(plyr)
#install.packages('xlsx')
library(openxlsx)
#install.packages('qqman')
#library(qqman)
#install.packages('ggplot2')
library(ggplot2)

library(readr)



#====================================

#####################################
#1. Specify Parameter settings
#####################################

#====================================
#===
#quality control settings
#===
mind<-1  #Mind: tolerated call rate of individuals => a mind of 0.10 corresponds to minimum 90% call-rate

#No LD pruning in this script (as recommended by Meyermans et al. 2020)

#===
#Specify ROH settings: divided in 2 classes (more info: https://www.cog-genomics.org/plink/1.9/ibd#homozyg)
#===

#=
#1. Window settings:
#=

#these settings are used to define window length to scan genome and
#assign 'scores' to individual snps. These scores determine if a snp is called in a 
#segment or not via threshold parameter

het<-0 #--homozyg-window-het AND --homozyg-het (it is important to specify both in PLINK; in this script het is used equal for both flags)
mis<-1 #--homozyg-window-missing
windowsnp<-0 ##--homozyg-window-snp => this parameter will be adjusted within script to parameter l!
windowthreshold<-0.05 #--homozyg-window-threshold => this parameter will be adjusted within script (t=floor(Nout+1L,3)!
Nout<-2 #The number of outer SNPs you don't want to allow for ROH-analysis (calculation of windowthreshold)

#=
#2.Segment settings:
#=

#After segment calling via window settings: extra requirements for segments!
kb<-1000 #--homozyg-kb: the minimum segment length in kb
l<-0    #--homozyg-snp => this parameter will be adjusted within script!
gap<-500 #--homozyg-gap: the maximal gap in kb between two SNPs within an ROH
density<- 60  #--homozyg-density #the minimal average density of an ROH (expressed as 1SNP/xxx kb, so in our example: minimal 1SNP/150kb)
####

####

#==
#Specify minimum number of animals that have to be present in population to perform ROH-analysis for that population (set to 1 if you want to do it for all populations)
#==

minimum_number_animals_population<-1


#====================================

#####################################
#2. Specify species to use and dataset
#####################################

#====================================



#===  
#  dataset name 
#===  

args <- commandArgs(trailingOnly = TRUE)

name <- args[1]
#name = "NE100_100cm_rep1"
if (grepl("_50", name)) {
data_name  <- sub("_50", "", name)
}else{
  data_name = name
}
#=====
#  Create directorie 
# =======

file.dir <- paste0("./Dataset/",data_name)

#=
#1. create output directory and directory to write  summary results of all breeds within species
#=

dir.create("./output")
dir.create(paste0("./output/",data_name))
dir.create("./output/all_chr_ROH")
dir.create(paste0("./output/all_chr_ROH/",data_name))

#=
#2. create map 'species' within ./output directory and ./output/all_chr_ROH directory
#=
if (grepl("_50", name)) {
dir.create(paste0("./output/",data_name,"/Half_Dens"))
dir.create(paste0("./output/all_chr_ROH/",data_name,"/Half_Dens"))

directory <- paste0("./output/",data_name,"/Half_Dens")

dir.create(paste0("./Dataset/",data_name,"/PLINK_out"))
dir.create(paste0("./Dataset/",data_name,"/PLINK_out/Half_Dens"))
dir.create(paste0("./Dataset/",data_name,"/PLINK_out/Half_Dens/Meyerman"))

plink.out <- paste0("./Dataset/",data_name,"/PLINK_out/Half_Dens/Meyerman")  

}else{
  #dir.create(paste0("./output/",data_name,"/Full_Dens"))
  #dir.create(paste0("./output/all_chr_ROH/",data_name,"/Full_Dens"))
  
  #directory <- paste0("./output/",data_name,"/Full_Dens")
  dir.create(paste0("./Dataset/",data_name,"/PLINK_out"))
  dir.create(paste0("./Dataset/",data_name,"/PLINK_out/Full_Dens"))
  dir.create(paste0("./Dataset/",data_name,"/PLINK_out/Full_Dens/Meyerman"))
  
  plink.out <- paste0("./Dataset/",data_name,"/PLINK_out/Full_Dens/Meyerman")   
}

#remove possible summary information with Froh from previous run
rm(summary_ROH_breeds)


#=====================
#If you did everything right, from this point on, you should not have to change any code anymore. 
#you can now just run the script
#=====================


#==
#Load fam file
#==

fam<-read.table(paste0(file.dir,"/","beforeQC.fam"),h=F)
fam<-fam[1:2]#only keep first two columns
colnames(fam)<-c("FID","IID")#rename columns
# 



#===
#Load map file, determine amount of chromosomes (no sex chromosomes!), and length of genome (in kb) covered by SNPs (sum of distance between first and last SNP per chromosome)
#===

#load bim file: determine number of chromosomes
if (grepl("_50", name)) {
  map <- read.table(paste0(file.dir,"/beforeQC_50.bim"),h=F)
}else{
  map <-read.table(paste0(file.dir,"/beforeQC.bim"),h=F)
}
 
colnames(map)<-c("CHR","SNP_name","POS","BP","A1","A2")    #adapt colnames
map_chr<-max(map$CHR)    #Determine number of chromosomes without sex chromosomes

#Calculate length in bp per chromosome and length of total genome covered (in kb)
genome_length = 0
#for loop per chromosome
for (ks in 1:map_chr){
  #map file of chromosome i
  tmp_map<-map[map$CHR==ks,]
  #determine genome length of this chromosome and add to total genome length
  genome_length<-genome_length+((max(tmp_map$BP)-min(tmp_map$BP))/1000) #divide by thousand to get kb
}



# #start loop for each population
# # for (k in 1:nrow(pop)){
# #   
#   #Print the current population
#   print(paste(species,dataset,pop[k,1]),sep=" ")
#   
#   #At beginning of run: remove files in directory (from possible previous runs)! 
#   #=>if error in next run, do not proceed with previous
#   unlink(x="./Data/ROH_analyse_population*", recursive = FALSE)
#   unlink(x="./Data/DummyGenotypes_ROH*", recursive = FALSE)
#   
#   #Choose population to use for ROH-analysis:
#   #Give this population a standard name and format: ROH_analyse_x
#   chosen_pop<-as.data.frame(as.character(pop[k,1]))
#   write.table(chosen_pop,paste("./Data/",species,"/Keep_FID.txt",sep=""),row.names = F,col.names = F,quote=FALSE)

 
#Load bed and bim beforeQC

if (grepl("_50", name)) {
  system(paste0(plink," --chr-set ",map_chr," --bfile ",file.dir,"/beforeQC_50","  --allow-extra-chr --make-bed --out ",plink.out,"/beforeQC"))

  }else{
  system(paste0(plink," --chr-set ",map_chr," --bfile ",file.dir,"/beforeQC","  --allow-extra-chr --make-bed --out ",plink.out,"/beforeQC"))
  
}



#Load QC-controlled datafile and make ped and map
system(paste0(plink," --chr-set ",map_chr," --bfile ",plink.out,"/beforeQC --allow-extra-chr --make-bed --recode --out ",plink.out,"/ROH_analyse"))


# name of Chromosome 
chr <- 1


  #Print species, dataset and population and chromosomes that are under evaluation
  print(paste("-data_name-",data_name,"-chromosome-",chr,sep=" "))
  
  
  #======
  #1.2.2 Define dataset
  #======
  
  #===
  #define data frame to fill in during for loop: general statistics populational level
  #===
  
  #Create empty matrix to later fill in with results
  Froh<-data.frame(matrix(ncol=18,nrow=1))
  #Colnames
  names(Froh)<-c("run","Het","Nanimals","NROH","Nsnps","KBROH","dummy_length","FROH","FROH1_2","FROH2_4","FROH4_8","FROH8_16","FROH16_","FROH5_","Fgrm","Fhom","Funi","l")
  
  #Write  details of number allowed heterozygotes and missing snps in ROH
  Froh[1,c("run")]<-paste(het,"_het_",mis,"_mis",sep="")
  
  #Write  number of animals in population
  Froh[1,c("Nanimals")]<-1
  
  #Load specific fam file for chosen population
  fam<-read.table(paste0(plink.out,"/ROH_analyse.fam"),h=F)
  colnames(fam)<-c("FID","IID")
  
  
  #===
  #define data frame to fill in during for loop: per animal statistics
  #===
  
  #Make empty data frame to fill in:
  ROH_animal<-data.frame(matrix(ncol=2,nrow=nrow(fam)))
  names(ROH_animal)<-c("FID","IID")
  
  #FID and IID information
  ROH_animal[,1:2]<-fam[,1:2]
  
  
  #======
  #1.2.3. ROH-analyse on whole dataset without pruning:
  #======
  
  #===
  #1. calculate Mean heterozygosity population using --hardy in plink
  #===
  
  system(paste0(plink," --chr-set ",map_chr," --allow-extra-chr --bfile ",plink.out,"/ROH_analyse --chr ",chr," --hardy --out ",plink.out,"/ROH_analyse"))
  
  #Load results
  plink.hwe <- read.table(file=paste0(plink.out,"/ROH_analyse.hwe"),fill=TRUE,header=TRUE)
  #write  observed heterozygosity
  heter<-mean(plink.hwe$O.HET.,na.rm = T) #
  
  #Write  mean heterozygosity value in data frame
  Froh[c(1),c("Het")]<-mean(plink.hwe$O.HET.,na.rm = T)
  
  #===
  #2. Mean heterozygosity per animal: --het in plink
  #===
  
  system(paste0(plink," --chr-set ",map_chr," --allow-extra-chr --bfile ",plink.out,"/ROH_analyse --chr ",chr," --het --out ",plink.out,"/ROH_analyse"))
  
  #load results
  plink.het <- read.table(file=paste0(plink.out,"/ROH_analyse.het"),fill=TRUE,header=TRUE)
  #Calculate observed heterozygosity per animal
  plink.het$O.Het<-(plink.het$N.NM.-plink.het$O.HOM.)/plink.het$N.NM.
  
  #Write  results
  ROH_animal$O_het<-plink.het$O.Het
  
  #===
  # 3 Number of SNPs per window
  #===
  
  #l=log(0.05/(number_of_SNPs*number of animals))/log(1-mean_heterozygosity)
  l <- (log(0.05/(nrow(map)*nrow(fam)))/log(1-heter)) #Lencz et al. 2007 Lencz, T., Lambert, C., DeRosse, P., Burdick, K. E., Morgan, T. V., Kane, J. M., ... & Malhotra, A. K. (2007). Runs of homozygosity reveal highly penetrant recessive loci in schizophrenia. Proceedings of the National Academy of Sciences, 104(50), 19942-19947.
  #or Purfield et al.2012 (Purfield, D. C., Berry, D. P., McParland, S., & Bradley, D. G. (2012). Runs of homozygosity and population history in cattle. BMC genetics, 13(1), 70.)
  #Round l parameter to an integer number
  l<-round(l)
  if (is.na(l)){l<-50} #In a rare occasion l is set to NA, set l = 50 in this case
  print(l) 
  
  #Change window-snp to l
  windowsnp<-l ##--homozyg-window-snp
  
  #Write  results
  Froh[c(1),c("l")]<-l
  
  #Calculate windowthreshold based on formula from Meyermans&Gorssen et al, 2020:
  windowthreshold<-floor(1000*((Nout+1)/l))/1000
  print(windowthreshold) 
  
  
  
  #===
  #4. ROH analysis in Plink (--homozyg)
  #===
  
  #ROH command in PLINK with all parameters as specified at start of script
  system(paste0(plink," --chr-set ",map_chr," --allow-extra-chr --bfile ",plink.out,"/ROH_analyse --chr ",chr," --homozyg group --homozyg-window-het ",het," --homozyg-het ",het," --homozyg-window-missing ",mis," --homozyg-snp ",l," --homozyg-kb ",kb," --homozyg-window-snp ",windowsnp,"  --homozyg-window-threshold ",windowthreshold,"  --homozyg-density ",density," --homozyg-gap ",gap," --out ",plink.out,"/ROH_analyse"))
  
  command<- paste0(plink," --chr-set ",map_chr," --allow-extra-chr --bfile ",plink.out,"/ROH_analyse --chr ",chr," --homozyg group --homozyg-window-het ",het," --homozyg-het ",het," --homozyg-window-missing ",mis," --homozyg-snp ",l," --homozyg-kb ",kb," --homozyg-window-snp ",windowsnp,"  --homozyg-window-threshold ",windowthreshold,"  --homozyg-density ",density," --homozyg-gap ",gap," --out ",plink.out,"/ROH_analyse")
  
  readr::write_lines(command, file = paste0(plink.out,"/PLINK_mey_commands.txt") )
  

  #================
  ## Inbreeding as Fhat1: --ibc command in plink (~variance-standardized relationship minus 1, Fhat2 (~Plink) and Fhat3 (~Yang) 
  #================
  
  #Run command for different inbreeding values via PLINK
  system(paste(plink," --chr-set ",map_chr," --allow-extra-chr --bfile ",plink.out,"/ROH_analyse --chr ",chr,"  --out ",plink.out,"/ROH_analyse",sep=""))
  

  
















