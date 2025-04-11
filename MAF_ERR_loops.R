# =============
# PLINK in R studio
# =============

source("/usr/share/lmod/lmod/init/R")
Sys.setenv(MODULEPATH = '/cluster/modules/all')
module("load PLINK/1.9b_6.17-x86_64")


# =======================
# Libraries 
.libPaths("/mnt/users/odwa/R/myLib" )


library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(readr)
library(stringi)
library(doParallel)
library(foreach)

source("~/paper-1/ROH-vs-IBD/v2_ROH_detection_functions.R")
source("~/paper-1/ROH-vs-IBD/vcferr_spTallman.R")
source("/mnt/users/odwa/paper-1/ROH-vs-IBD/F_in_ROH_functions.R")



# ==============
#     Cluster
# =============
num_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

# Use one less than available cores
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clusterEvalQ(cl, {
  .libPaths("/mnt/users/odwa/R/myLib" )
  library(data.table)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringi)
  
  source("/usr/share/lmod/lmod/init/R")
  Sys.setenv(MODULEPATH = '/cluster/modules/all')
  module("load PLINK/1.9b_6.17-x86_64")
  
  source("~/paper-1/ROH-vs-IBD/IBD_vs_ROH_functions.R")
  source("~/paper-1/ROH-vs-IBD/v2_ROH_detection_functions.R")
  source("~/paper-1/ROH-vs-IBD/vcferr_spTallman.R")
  source("/mnt/users/odwa/paper-1/ROH-vs-IBD/F_in_ROH_functions.R")
  
})

# =============
# Path parameters and combo lists
# =============

NEs = c("NE100","NE200")
Reps = c(1:10)
Densities = c("Full_Dens","Half_Dens")
Methods <- c("Meyermans", "Default", "Norm", "Norm_small")
Errors = c(0,0.005,0.01,0.02)
MAFs = c(0 , 0.01,0.05)
Gens = c(10,20,50,100)
proj.dir = "/mnt/project/SimData/Paper_1"

combinations <-  expand.grid(NE = NEs,Rep = Reps, Dens = Densities, Method = Methods, Err = Errors, MAF = MAFs, Gen = Gens)

# "
# MAF on-top of of error
# "

foreach(i = 1:nrow(combinations)) %dopar% {
  NE = combinations$NE[i]
  Rep = combinations$Rep[i]
  Dens = combinations$Dens[i]
  Method = combinations$Method[i]
  Err = combinations$Err[i]
  MAF = combinations$MAF[i]
  Gen = combinations$Gen[i]

  data.dir<- paste0(proj.dir,"/",NE,"/Rep",Rep)
  dens.dir <-paste0(data.dir,"/",Dens)
  err.dir <- paste0(dens.dir,"/MAF_and_Error")

  out.dir = paste0(err.dir,"/",Gen,"_",Err,"_",MAF)
  dir.create(out.dir, recursive = TRUE)



  fread(paste0(data.dir,"/pedigree_rel_generations.txt")) %>%
    filter(birth == Gen) %>%
    mutate(FID = 0) %>%
    select(FID,id) %>%
    write.table(file = paste0(out.dir,"/",Gen,"_rel_inds.txt"),
                row.names = F,
                col.names = F,
                quote = F,
                sep = '\t')

  system(paste0("plink --bfile ",dens.dir,"/beforeQC",
                " --keep ",out.dir,"/",Gen,"_rel_inds.txt  --recode vcf --const-fid --out ",out.dir,"/",Gen,"_inped_1"))
  # Run the AWK command within R using system()
  system(paste0(
    "awk 'BEGIN { OFS=\"\\t\" } { for (i = 10; i <= NF; i++) { gsub(\"/\", \"|\", $i) } print }' ",
    out.dir,"/",Gen,"_inped_1.vcf", " > ",out.dir,"/",Gen,"_inped_2.vcf"
  ))

  # Generate the errors
  vcferr(paste0(out.dir,"/",Gen,"_inped_2.vcf"), Err, paste0(out.dir,"/",Gen,"_")  )

  system(paste0("grep '\\#' ",out.dir,"/",Gen,"_inped_2.vcf | cat - ",out.dir,"/",Gen,"_",Err,".err.vcf > ",out.dir,"/",Gen,"_",Err,".h.err.vcf"))


  if (MAF == 0) {
    system(paste0("plink --vcf ",out.dir,"/",Gen,"_",Err,".h.err.vcf",
                  " --make-bed --out ", out.dir,"/inped_3"))
  }else{
  system(paste0("plink --vcf ",out.dir,"/",Gen,"_",Err,".h.err.vcf",
                " --maf ",MAF," --make-bed --out ", out.dir,"/inped_3"))
  }
  inped = paste0("inped_3")

  # Detect the ROH
  Gorssen(input.dir = out.dir, inped =inped ,plink.out = paste0(out.dir,"/Meyermans"))
  Default(input.dir = out.dir,inped = inped ,plink.out = paste0(out.dir,"/Default"))
  Norm(input.dir = out.dir,inped = inped ,plink.out = paste0(out.dir,"/Norm"))
  Norm_small(input.dir = out.dir,inped = inped ,plink.out = paste0(out.dir,"/Norm_small"))

} # dataset



# ===================
#
#   Benchmarking
#
# ===================

# Generate all combinations of dataset, density, method, and err
combinations <- expand.grid(NE= NEs ,Rep = Reps, Dens = Densities, Method = Methods, Err = Errors, MAF = MAFs)

# Parallel processing over all combinations


# Parallel loop
foreach(i = 1:nrow(combinations)) %dopar% {

  # Get the current combination
  NE <- combinations$NE[i]
  Rep <- combinations$Rep[i]
  Dens <- combinations$Dens[i]
  Method <- combinations$Method[i]
  Err <- combinations$Err[i]
  MAF <- combinations$MAF[i]

  data.input.dir <- paste0(proj.dir, "/", NE,"/Rep",Rep)

  marker_pos <- read_csv(paste0(data.input.dir, "/markerPosition_kb.txt"), show_col_types = FALSE)
  true_IBD <- read_csv(paste0(data.input.dir, "/All_Ind_true_IBD.res.bz2"), show_col_types = FALSE)
  IBD_pos <- read.csv(paste0(data.input.dir, "/true_IBD_position.txt"))
  ped <- read_table(paste0(data.input.dir, "/pedigree_rel_generations.txt"))
  obs_IBD <- read_csv(paste(data.input.dir, "/extend_bin_ADAM_ROH.res.bz2", sep = ""))
  obs_homo <- read_csv(paste0(data.input.dir, "/all_ind_observed_homozyg_QTL.res.bz2"))

  info_ADAM_ROH_bin <- read.csv(paste0(data.input.dir, "/info_ADAM_roh_bin.txt")) %>%
    filter(ID %in% ped$id) %>%
    mutate(start_kb = start_bpp / 1000,
           stop_kb = stop_bpp / 1000)

  ROH_inf <- data.frame()
  for (Gen in Gens) {

    roh.dir <- paste0(data.input.dir,"/",Dens,"/MAF_and_Error/",Gen,"_",Err,"_",MAF,"/",Method)

    ROH_inf_t <- read_table(paste0(roh.dir, "/ROH_analyse.hom")) %>%
      select(c(IID, SNP1, SNP2, POS1, POS2, KB, NSNP)) %>%
      rename(ID = IID) %>%
      filter(ID %in% ped$id) %>%
      group_by(ID) %>%
      mutate(run = row_number())
    if (nrow(ROH_inf_t) > 0) {
      ROH_inf <- rbind(ROH_inf, ROH_inf_t)
    }

  }

  PLINK_ROH_bin <- RohBinGeno(ROH_inf, marker_pos)

  RohPower(info_ADAM_ROH_bin, PLINK_ROH_bin, marker_pos, tmp.dir = roh.dir)
  RohTruePos(ROH_inf, obs_IBD, obs_homo, marker_pos, tmp.dir = roh.dir)

} # End of foreach loop

# Stop the cluster after processing


#

# Generate all combinations of dataset, density, method, and err
combinations <-  expand.grid(NE = NEs,Rep = Reps, Dens = Densities, Method = Methods, MAF = MAFs, Gen = Gens)

results_list <- foreach(i = c(1:nrow(combinations)), .packages = c("foreach", "doSNOW","tidyr","data.table"),
                        .export = c("process_SNP","hom_in_ROH","F_in_ROH")) %dopar% {
                          NE = combinations$NE[i]
                          Rep = combinations$Rep[i]
                          Dens = combinations$Dens[i]
                          Method = combinations$Method[i]
                          Gen = combinations$Gen[i]
                          MAF = combinations$MAF[i]
                          
                          data.dir = paste0(proj.dir,"/",NE,"/Rep",Rep)
                          roh.dir <- paste0(data.dir,"/",Dens,"/MAF_and_Error/",Gen,"_0_",MAF,"/",Method)
                          

                          ped <- fread(file.path(data.dir, "pedigree_rel_generations.txt")) %>%
                            filter(birth == Gen)

                          roh <- fread(paste0(roh.dir,"/ROH_analyse.hom")) %>%
                            filter(IID %in% ped$id)

                          marker_freq_all <- read.csv(paste0(data.dir,"/AlleleFreqs_gens.txt" )) %>%
                            filter(GEN == Gen)
                          homozyg_marker_all <- read_csv(paste0(data.dir,"/",Dens,"/homozyg_test_markers.txt"),show_col_types = FALSE) %>%
                            filter(Generation == Gen)

                          if (nrow(roh)>0) {
                            hom_in_ROH(data.dir, roh.dir, Dens, Gen)

                            t2 <- F_in_ROH(homozyg_marker_all, marker_freq_all, data.dir, roh.dir,Gen)
                            t2 <- t2 %>%
                              mutate(NE = NE,
                                     Rep = Rep,
                                     Dens = Dens,
                                     Method = Method,
                                     Gen = Gen,
                                     MAF = MAF)

                            return(t2)
                          }

                        }

# Stop the cluster after processing

stopCluster(cl)

ID_result <- do.call(rbind, results_list)


write.table(ID_result, file= file.path(proj.dir,"Results/MAF_F_in_ROH_results"),sep = "," , quote = FALSE, row.names = FALSE)


