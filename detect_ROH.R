

# ===============
#   Making PLINK runable from Rstudio 
# ===============

source("/usr/share/lmod/lmod/init/R")
Sys.setenv(MODULEPATH = '/cluster/modules/all')
module("load PLINK/1.9b_6.17-x86_64")

source("~/paper-1/ROH-vs-IBD/v2_ROH_detection_functions.R")


#NEs = c("NE200")
NEs = c("NE100","NE200")

Reps = c(1:10)
Densities = c("Full_Dens", "Half_Dens")

proj.dir <- paste0("/mnt/project/SimData/Paper_1")



for (NE in NEs) {
  for (Rep in Reps) {
    for (Dens in Densities) {
      
    plink.path = paste0(proj.dir,"/",NE,"/Rep",Rep,"/",Dens)
    #function(input.dir, inped,plink.out)
    # Detect the ROH
    Gorssen(input.dir = plink.path,inped = "beforeQC" ,plink.out = paste0(plink.path,"/Meyermans"))
    Default(input.dir = plink.path,inped = "beforeQC" ,plink.out = paste0(plink.path,"/Default"))
    Norm(input.dir = plink.path,inped = "beforeQC" ,plink.out = paste0(plink.path,"/Norm"))
    Norm_small(input.dir = plink.path,inped = "beforeQC" ,plink.out = paste0(plink.path,"/Norm_small"))
    
    }
  }
}

