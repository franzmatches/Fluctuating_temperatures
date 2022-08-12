#########################################################################################
### DISPERSAL TEST LONG VS SHORT CORRIDORS #####
#########################################################################################

#preable
library(openxlsx)
library(tidyverse) # data wrangling and plotting
library(ggplot2)


#load data
data<-read.xlsx("Data/Dispersal_experiments_data.xlsx", sheet = "Prey_species", colNames = T)

#perform Mann-Whitney-Wilcoxon Test
wilcox.test(N_individuals_2h ~ Corridor,data = data %>% filter(Species == "P.caudatum"))
wilcox.test(N_individuals_2h ~ Corridor,data = data %>% filter(Species == "B.japonicum"))
wilcox.test(N_individuals_2h ~ Corridor,data = data %>% filter(Species == "S.teres"))
wilcox.test(N_individuals_2h ~ Corridor,data = data %>% filter(Species == "Colpidium"))

#plot boxplots to see difference between long and short facet species
ggplot(data, aes (x = Corridor, y = N_individuals_2h, color = Species))+
  geom_boxplot()+
  facet_wrap(~Species, scales = "free")+
  theme_bw()+
  ggtitle("Prey species")
ggsave("dispersal_boxplots.pdf", units = "in", width = 10, height = 10)
