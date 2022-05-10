### Figure 3 ###

## Preamble
library(ggh4x)
library(tidyverse) # data wrangling and plotting

source("Code/data_preparation.R")
palette <- c("#440154FF", "#3B528BFF","#21908CFF", "#5DC863FF")

#create a sub-dataset with just the last day of experiment
data_pooled_lastday <- data_pooled %>%
  filter(NumDays ==max(NumDays))
#select just the last day to analyse the Beta diversity at the end of the experiment
data_patches_lastday <- data_patches %>% 
  filter(NumDays ==max(NumDays) & Patch == "A")


## Shannon Diversity at end of experiment ##
shannon_end <- ggplot(data_pooled_lastday, aes( x = Temp_Regime, y = Sh_div, col = Temp_Regime))+
  geom_boxplot()+
  xlab("Temperature regime") +
  ylab("Shannon diversity") +
  theme_classic() +
  theme(legend.position = "none",
        aspect.ratio = 1,
        axis.text.x = element_text(size = 8)) +
  scale_colour_manual(values= palette) +
  scale_x_discrete(labels = c("Constant", "Fluctuating \n asynchronous",
                              "Fluctuating \n synchronous", "Static difference"))


## Beta Diversity between the two patches at end of experiment ##
beta_end <- ggplot(data_patches_lastday, aes( x = Temp_Regime, y = beta, col = Temp_Regime))+
  geom_boxplot()+
  xlab("Temperature regime") +
  ylab("Beta diversity") +
  #  labs(colour = "Temperature regime") +  
  theme_classic()+
  theme(legend.position = "none",
        aspect.ratio = 1,
        axis.text.x = element_text(size = 8)) +
  scale_colour_manual(values= palette) +
  scale_x_discrete(labels = c("Constant", "Fluctuating \n asynchronous",
                              "Fluctuating \n synchronous", "Static difference"))  


## Combine figures ##

end_plots <- ggpubr::ggarrange(shannon_end, beta_end,  labels = "auto", label.x = 0.08, label.y = 0.71)
ggsave("Results/Figure3.tiff", end_plots,units = "in", width = 10, height = 10)
