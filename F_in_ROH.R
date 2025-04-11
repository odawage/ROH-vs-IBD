
# source("/usr/share/lmod/lmod/init/R")
# Sys.setenv(MODULEPATH = '/cluster/modules/all')
# module("load PLINK/1.9b_6.17-x86_64")

# =============
#   libs 
# =============
.libPaths("/mnt/users/odwa/R/myLib")
library(foreach)
library(doParallel)
library(doSNOW)
library(data.table)
library(readr)
library(tidyverse)

source("/mnt/users/odwa/paper-1/ROH-vs-IBD/F_in_ROH_functions.R")

# =============
#   Pathing parameters 
# =============


NEs = c("NE100","NE200")
Reps = c(1:10)
Densities <- c("Full_Dens", "Half_Dens")
Methods <- c("Meyermans", "Default", "Norm", "Norm_small")
Gens <- c(10,20,50,100)

proj.dir = "/mnt/project/SimData/Paper_1"


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
  library(readr)
  library(tidyverse)
})

# ==============
# calculate allele freqs
# =============
df_list <- expand.grid(NE = NEs,Rep = Reps)

foreach(i = 1:nrow(df_list)) %dopar% {
  NE = df_list$NE[i]
  Rep = df_list$Rep[i]
  data.dir = paste0(proj.dir,"/",NE,"/Rep",Rep)
  
  Allele_freq_gen(data.dir, Gens)

  }

# ==============
# calculate frac homozyg
# =============
dens_list <- expand.grid(NE = NEs,Rep = Reps, Dens = Densities )

foreach(i = 1:nrow(dens_list)) %dopar% {
  NE = dens_list$NE[i]
  Rep = dens_list$Rep[i]
  Dens = dens_list$Dens[i]
  data.dir = paste0(proj.dir,"/",NE,"/Rep",Rep)
  
  homo_test(data.dir, NE, Rep, Dens, rel_gen)
}



# Generate all combinations of dataset, density, method, and err
combinations <- expand.grid(NE = NEs,Rep = Reps, Dens = Densities, Method = Methods, Gen = Gens)

results_list <- foreach(i = c(1:nrow(combinations)), .packages = c("foreach", "doSNOW","tidyverse","data.table"),
                        .export = c("process_SNP","hom_in_ROH","F_in_ROH")) %dopar% { 
  
  NE = combinations$NE[i]
  Rep = combinations$Rep[i] 
  Dens = combinations$Dens[i]
  Method = combinations$Method[i]
  Gen = combinations$Gen[i]
  
  data.dir = paste0(proj.dir,"/",NE,"/Rep",Rep)
  roh.dir <-  paste0(data.dir,"/",Dens,"/",Method)

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
             Gen = Gen)
    
    return(t2)
  }
  
}

# Stop the cluster after processing

stopCluster(cl)

ID_result <- do.call(rbind, results_list)


write.table(ID_result, file= file.path(proj.dir,"Results/F_in_ROH_results"),sep = "," , quote = FALSE, row.names = FALSE)

