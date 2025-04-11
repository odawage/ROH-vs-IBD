

.libPaths("/mnt/users/odwa/R/myLib")

library(readr)
library(data.table)
library(tidyverse)

# 
source("/usr/share/lmod/lmod/init/R")
Sys.setenv(MODULEPATH = '/cluster/modules/all')
module("load PLINK/1.9b_6.17-x86_64")

source("~/paper-1/ROH-vs-IBD/v2_ROH_detection_functions.R")



NEs = c("NE100","NE200")
Reps = c(1:10)
Densities <- c("Full_Dens", "Half_Dens")
Methods <- c("Meyermans", "Default", "Norm", "Norm_small")

proj.dir = "/mnt/project/SimData/Paper_1"

comb_list <- expand.grid(NE = NEs,Rep = Reps, Dens = Densities)

for (i in 1:nrow(comb_list)) {
  NE = comb_list$NE[i]
  Rep = comb_list$Rep[i]
  Dens = comb_list$Dens[i]
  
  data.dir = paste0(proj.dir,"/",NE,"/Rep",Rep,"/",Dens)
  out.dir = file.path(data.dir,"FROH_expanded")
  dir.create(out.dir,recursive = TRUE)
  
  if (Dens == "Half_Dens") {
    markers <-fread(paste0(data.dir,"/MarkerQtlGenotypesRep1.res.bz2"))
    info <- data.frame(matrix(nrow = nrow(markers), ncol = 6, 0))
    info[ , 2] <- as.numeric(markers$ID)
    info <- unite(info, info, c(1:6), sep = " ")
    info$Geno <- markers$Geno
    write.table(info, file= paste0(data.dir,"/full.ped"), row.names = FALSE, col.names = FALSE,sep = " ",quote = FALSE )
  }else{
    file.copy(paste0(proj.dir,"/",NE,"/Rep",Rep,"/full.ped"), paste0(data.dir), overwrite = TRUE )
  }
  
  system(paste0("plink --file ",data.dir,"/full --make-bed --out ",out.dir,"/beforeQC"))    
  Gorssen(input.dir = out.dir,inped = "beforeQC" ,plink.out = paste0(out.dir,"/Meyermans"))
  Default(input.dir = out.dir,inped = "beforeQC" ,plink.out = paste0(out.dir,"/Default"))
  Norm(input.dir = out.dir,inped = "beforeQC" ,plink.out = paste0(out.dir,"/Norm"))
  Norm_small(input.dir = out.dir,inped = "beforeQC" ,plink.out = paste0(out.dir,"/Norm_small"))
  
}

# ================
#
# ================
genome_l <- data.frame()
for (i in 1:nrow(comb_list)) {
  NE = comb_list$NE[i]
  Rep = comb_list$Rep[i]
  Dens = comb_list$Dens[i]
  
  data.dir = paste0(proj.dir,"/",NE,"/Rep",Rep,"/",Dens)
  map <- fread(file.path(data.dir,"full.map"))
  df <- data.frame(first = min(map$V4), last = max(map$V4), length_bpp =(max(map$V4)-min(map$V4)),  NE = NE, Rep = Rep, Dens = Dens)
  genome_l <- rbind(genome_l, df)
  }



# ================
## calculate FROH
# ================
comb_list <- expand.grid(NE = NEs,Rep = Reps, Dens = Densities, Method = Methods )
ind_ROH <- data.frame()
for (i in 1:nrow(comb_list)) {
  NE = comb_list$NE[i]
  Rep = comb_list$Rep[i]
  Dens = comb_list$Dens[i]
  Method = comb_list$Method[i]
  
  roh.dir <- paste0(proj.dir,"/",NE,"/Rep",Rep,"/",Dens,"/FROH_expanded/",Method)

  
  ped <- fread(paste0(proj.dir,"/",NE,"/Rep",Rep,"/animalsRep1.res.bz2")) 
    
  roh <- fread(file.path(roh.dir,"ROH_analyse.hom.indiv")) %>% 
    dplyr::mutate(NE = NE, 
                  Rep = Rep,
                  Dens = Dens, 
                  Method = Method) %>% 
    dplyr::rename(ID = "IID")
  roh <- left_join(roh, ped[,c("id", "birth")], by = c("ID" = "id"))
  
  ind_ROH <- rbind(ind_ROH, roh)
  }
  
# ================
#
# ================
ind_IBD <- data.frame()
comb_list <- expand.grid(NE = NEs,Rep = Reps)
for (i in 1:nrow(comb_list)) {
  NE = comb_list$NE[i]
  Rep = comb_list$Rep[i]

  data.dir = paste0(proj.dir,"/",NE,"/Rep",Rep)
  
  IBD_seg <- fread(file.path(data.dir,"info_ADAM_roh_bin.txt"))
  
  IBD_seg <- IBD_seg %>% 
    mutate(length_bpp = stop_bpp-start_bpp) %>% 
    group_by(ID) %>% 
    dplyr::summarise(NSEG = n(),
                      KB = sum(length_bpp/1000),
                      KBAVG = mean(length_bpp/1000)
                      ) %>% 
    dplyr::mutate(NE = NE, 
                  Rep = Rep,
                  Dens = "Full_Dens", 
                  Method = "IBD")
  ped <- fread("/mnt/project/SimData/Paper_1/NE100/Rep1/animalsRep1.res.bz2")
  IBD_seg <- left_join(IBD_seg,  ped[,c("id", "birth")], by = c("ID"= "id"))
  
  ind_IBD <- rbind(ind_IBD, IBD_seg)
}

FROH_FIBD <- left_join(ind_ROH[,-c("FID","PHE")], ind_IBD[,1:6], by =c("NE", "Rep", "ID"), suffix = c("", "_IBD"))

FROH_FIBD <- left_join(FROH_FIBD, genome_l, by =c("NE", "Rep","Dens") )



steg_1 <- FROH_FIBD %>% 
  replace(is.na(.), 0) %>% 
  group_by( Method, NE, Dens, birth) %>%
  dplyr::summarise(m_FROH = mean(KB/(length_bpp/1000)),
                   m_FIBD = mean(KB_IBD/(length_bpp/1000)))
                   

steg_2 <- steg_1 %>% 
  group_by(Method, NE,Dens) %>% 
  mutate(d_FROH = (m_FROH - lag(m_FROH))/(1-lag(m_FROH)),
         d_FIBD = (m_FIBD - lag(m_FIBD))/(1-lag(m_FIBD))
  )


write.table(steg_2, file.path(proj.dir,"Results/Delta_F"), row.names = FALSE, col.names = TRUE,sep = " ",quote = FALSE)

# df <- steg_2 %>% 
#   filter(m_FROH != 0) %>% 
#   group_by(NE, Method, Dens) %>%
#   mutate(xpos = 25,
#          ypos =0.95 + ifelse(Dens =="Full_Dens",0.09,0 ))
# library(broom)
# 
# # Calculate regression equations and R-squared values for each group
# regressions <- df %>%
#   do(tidy(lm(-log(1-m_FROH) ~ birth, data = .))) %>%
#   select(!std.error:p.value) %>% 
#   pivot_wider(names_from = term, values_from = estimate) %>% 
#   mutate(eq = paste0("y = ",round(`(Intercept)`,2),"+x", round(`birth`,4) ),
#          beta = paste0("ΔF = ", round(`birth`,5) ),
#          ypos = 0.65 + ifelse(Dens =="Full_Dens",0.1,0 ))
# #filter(term != "(Intercept)") %>%
# 
# plot <- ggplot(df, aes(x=birth, y=-log(1-m_FROH), color=Dens)) +
#   geom_point(size = 0.1, alpha= 0.3) +
#   #geom_smooth(method=lm, se=FALSE, formula = y ~ x) +
#   geom_line()+
#   geom_line(aes(y = -log(1-m_FIBD)), colour = "black")+
#   facet_grid(NE ~ Method)+
#   geom_text(data = regressions, aes(x = 35, y = ypos, label = beta, color = Dens), 
#             size = 3) +
#   #scale_colour_manual(values = pal[c(5,3)])+
#   theme_minimal()+
#   theme( legend.position = "top") +
#   guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
#   labs(#title = "F(IBD) vs F(ROH) split by generation and setting",
#     y = expression(paste(-log(1-F[t]) )),
#     x = "Generation" , 
#     color = "Density")
# 
# plot
