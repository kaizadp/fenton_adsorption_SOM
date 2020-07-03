source("code/0-packages.R")

relabund = read.csv(GOETHITE_RELABUND)

relabund2 = 
  relabund %>% 
  filter(!adsorbed == "new") %>% 
  dplyr::select(Forest, fenton, Class, relabund) %>% 
  group_by(Forest, fenton, Class) %>% 
  dplyr::summarise(peaks = sum(relabund)) %>% 
  spread(Class, peaks)

relabund2$DV = as.matrix(relabund2[,3:8])

library(compositions)

man = manova(ilr(clo(DV)) ~ fenton, data = relabund2)
summary(man)

# -------
master = read.csv(FTICR_MASTER_LONG) 
master2 = 
  master %>% 
  spread(Mass, presence) %>% 
  replace(is.na(.),0)
master2$DV = as.matrix(master2[,3:3214])

man = manova(ilr(clo(DV)) ~ Treatment, data = master2)
summary(man)