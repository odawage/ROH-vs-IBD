
library(readr)
library(tidyverse)
library(pedigreemm)
library(data.table)

datasets = c("NE200_100cm_rep0", "NE200_100cm_rep1", "NE200_100cm_rep2", "NE100_100cm_rep1", "NE100_100cm_rep2")


ped_all <- data.frame()
for (dataset in datasets) {
  rep <- str_extract(dataset, "\\d+(?=$)")
  ped <- fread(paste0("/mnt/users/odwa/PLINK/Code/Dataset/",dataset,"/animalsRep",rep,".res.bz2") )%>%
                 mutate(across(-c(birth, death), ~ifelse(. == 0, NA, .)))
               
               ped <- pedigreemm::editPed(sire = ped$sire, dam = ped$dam, label = ped$id)
               
               ped_obj <- pedigreemm::pedigree(sire = ped$sire, dam = ped$dam, label = ped$label)
               
               ped$inb  <- pedigreemm::inbreeding(ped_obj)
               ped <- ped %>% mutate(Dataset = dataset,
                                     simNE = str_extract(dataset, "(NE\\d+)") )
               ped_all <- rbind(ped_all, ped)
}

library(broom)
# DF <- -log(1-Ft0), i use all individual inbreeding values in my 
# NeFt = 1/2ΔFt 
  
values <- ped_all %>% group_by(simNE) %>% 
  do(tidy(lm(-log(1-inb) ~ gene, data = .))) %>%
  select(!std.error:p.value) %>% 
  pivot_wider(names_from = term, values_from = estimate) %>% 
  mutate(NE = 1/(2*gene))


ped_all %>% ggplot(aes(x = gene, y= inb))+
  geom_point()+
  geom_abline(aes(data = values, intercept =`(Intercept)`, slope = gene) )+
  facet_wrap(~factor(simNE))

values <- values %>%
  mutate(label = glue("y = {round(`(Intercept)`, 5)} + {round(gene, 5)} * x\nNE = {round(NE, 5)}"))

ped_all %>%
  ggplot(aes(x = gene, y = inb)) +
  geom_point(size =0.3,aes(color = factor(simNE)) ) +
  geom_abline(
    data = values,
    mapping = aes(intercept = `(Intercept)`, slope = gene),
    size = 1
  )+
  geom_text(
    data = values,
    mapping = aes(x = -Inf, y = Inf, label = label),
    hjust = -0.1, vjust = 1.5, size = 4
  ) +
  facet_wrap(~factor(simNE)) +
  scale_color_discrete(name = "simNE") +
  theme_minimal() +
  labs(x = "Generation", y = "Inbreeding Coefficient")













