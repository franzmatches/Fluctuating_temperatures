### Supplementary Figures ###

## Preamble
library(ggh4x)
library(tidyverse) # data wrangling and plotting
library(dplyr)
library(plotrix) #for calculating standard error

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

#----------------------------------------------------------------------------------------
# Supplementary Figure 2 - alpha richness
#----------------------------------------------------------------------------------------

#group patches together to obtain average alpha richness 
data_average_alpha <- data_patches %>% ungroup () %>% 
  group_by(Replicate, Temp_Regime, Corridor, NumDays) %>% 
  summarize(average_alpha = mean(Sp_rich))

#generate average alpha richness and SE, put into a new data frame for plotting
data_alpha_plotting <- data_average_alpha %>% ungroup() %>%
  group_by(Temp_Regime, Corridor, NumDays) %>%
  dplyr :: summarise(mean_alpha = mean(average_alpha),
                     se_alpha = std.error(average_alpha))

alpha_plot <- 
  ggplot(data_alpha_plotting, aes(x = NumDays, y = mean_alpha, group = Corridor, colour = Corridor))+
  geom_point(position=position_dodge(0.2))+
  geom_line(position=position_dodge(0.2))+
  geom_errorbar(aes(ymin = mean_alpha - se_alpha,
                    ymax = mean_alpha + se_alpha),
                width=.5, size = .6,
                position = position_dodge(0.2))+
  facet_grid(Corridor ~ Temp_Regime, scales = "fixed",
             labeller = ggplot2::labeller(Temp_Regime = temp_labs))+
  xlab("Time (days)")+
  ylab("Average alpha richness")+
  scale_colour_manual(values = c("#440154FF","#21908CFF"))+
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black"))

ggsave("Results/Supplementary_Figure2.tiff", alpha_plot, units = "in", width = 10, height = 10)

#----------------------------------------------------------------------------------------
# Supplementary Figure 3 - gamma richness
#----------------------------------------------------------------------------------------

#generate average alpha richness and SE, put into a new data frame for plotting
data_gamma_plotting <- data_pooled %>% ungroup() %>%
  group_by(Temp_Regime, Corridor, NumDays) %>%
  dplyr :: summarise(mean_gamma = mean(Sp_rich),
                     se_gamma = std.error(Sp_rich))

gamma_plot <- 
  ggplot(data_gamma_plotting, aes(x = NumDays, y = mean_gamma, group = Corridor, colour = Corridor))+
  geom_point(position=position_dodge(0.2))+
  geom_line(position=position_dodge(0.2))+
  geom_errorbar(aes(ymin = mean_gamma - se_gamma,
                    ymax = mean_gamma + se_gamma),
                width=.5, size = .6,
                position = position_dodge(0.2))+
  facet_grid(Corridor ~ Temp_Regime, scales = "fixed",
             labeller = ggplot2::labeller(Temp_Regime = temp_labs))+
  xlab("Time (days)")+
  ylab("Gamma richness")+
  scale_colour_manual(values = c("#440154FF","#21908CFF"))+
  theme_classic()+
  theme(strip.background = element_rect(colour = "black", fill = "white", linetype = "blank"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        aspect.ratio = 1,
        legend.position = "none",
        panel.border = element_rect(fill = NA, colour = "black"))

ggsave("Results/Supplementary_Figure3.tiff", gamma_plot, units = "in", width = 10, height = 10)

