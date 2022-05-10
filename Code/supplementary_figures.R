### Supplementary Figures ###

## Preamble
library(ggh4x)
library(tidyverse) # data wrangling and plotting

source("Code/data_preparation.R")

#----------------------------------------------------------------------------------------
# Supplementary Figure 1
#----------------------------------------------------------------------------------------

#pivot for obtaining the species variable for plotting
data_pivot <- data_patches%>%
  pivot_longer(cols = -c(Date, Replicate, Temp_Regime,
                         Corridor,Patch, beta, NumDays),
               names_to ="species",
               values_to= "abundance")%>%
  group_by(Temp_Regime, Corridor,Patch,Replicate)


data_merged_corridors<-data_pivot %>% ungroup () %>% 
  group_by(Date, Replicate, Temp_Regime,Patch, NumDays, species) %>% 
  summarize(summed_abundance = sum(abundance))


temp_labs <- c("Constant", "Fluctuating\nasynchronous", "Fluctuating\nsynchronous", "Static difference")
names(temp_labs) <- c("Constant", "Fluctuating_Asynchro", "Fluctuating_Synchro", "Static_Diff")

abundance_plot <- ggplot(data_pivot, aes(x = NumDays, y = log(abundance+1), 
                                         colour = species, fill = species))+
  geom_smooth(method = "loess",alpha = 0.15)+
  facet_nested(Corridor + Patch ~ Temp_Regime,scales = "fixed",
               strip = ggh4x::strip_nested(size="constant",bleed=T),
               labeller = ggplot2::labeller(Temp_Regime = temp_labs))+
  xlab("Time (days)")+
  ylab("Log abundance")+
  labs(colour = "Species", fill = "Species") +
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"))+
  scale_colour_viridis_d(labels = c("Blepharisma", "Bursaria", "Colpidium", "Didinium", "Homalozoon", "Paramecium", "Spirostomum"))+
  scale_fill_viridis_d(labels = c("Blepharisma", "Bursaria", "Colpidium", "Didinium", "Homalozoon", "Paramecium", "Spirostomum"))

ggsave("Results/Supplementary_Figure1.tiff", abundance_plot, units = "in", width = 10, height = 10)
