### Figure 2 ###

## Preamble
require(ggeffects) # marginal effects extraction
require(ggplot2) # plotting
require(glmmTMB) # model fitting

source("Code/data_preparation.R")
palette <- c("#440154FF", "#3B528BFF","#21908CFF", "#5DC863FF") #colour paletter

### Fit the final species diversity model
shannon_mod_3 <- glmmTMB::glmmTMB(Sh_div~ splines::ns(NumDays,3) + Temp_Regime + Corridor + 
                                    splines::ns(NumDays,3):Temp_Regime + splines::ns(NumDays,3):Corridor +
                                    ar1(as.factor(NumDays) + 0 | Replicate) + (1|Replicate), 
                                  data = data_pooled,family = "gaussian",REML=F)


### Species diversity across temperature regimes model ###
gg.df_Cor <- ggeffects::ggpredict(shannon_mod_3, c("NumDays[all]","Corridor[all]"))%>%
  left_join(attr(ggeffects::ggpredict(shannon_mod_3, c("NumDays[all]","Corridor[all]")), "rawdata", exact = TRUE))

temperature_time <- ggplot(gg.df_TEMP,aes(x=x,y=predicted,colour = group)) +  
  geom_line()+
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high,group=group),colour = NA, fill = "grey20",alpha=0.1)+
  geom_point(aes(y=response),alpha = 1, pch = 1)+
  xlab("Time (days)") +
  ylab("Landscape diversity") +
  labs(colour = "Temperature regime", fill = "Temperature regime") +  
  theme_classic()+
  theme(legend.position = "top",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        aspect.ratio = 1) +
  scale_colour_manual(values= palette,
                      labels = c("Constant", "Fluctuating\nasynchronous",
                                 "Fluctuating\nsynchronous", "Static\ndifference")) 

### Species diversity across corridor lengths model ###
gg.df_TEMP <- ggeffects::ggpredict(shannon_mod_3, c("NumDays[all]", "Temp_Regime[all]")) %>%
  left_join(attr(ggeffects::ggpredict(shannon_mod_3, c("NumDays[all]", "Temp_Regime[all]")), "rawdata", exact = TRUE))

corridor_time <- ggplot(gg.df_Cor,aes(x=x,y=predicted,colour = group)) +  
  geom_line()+
  geom_ribbon(aes(ymin = conf.low,
                  ymax = conf.high,group=group),colour = NA, fill = "grey20",alpha=0.1)+
  geom_point(aes(y=response),alpha = 1, pch = 1)+
  xlab("Time (days)") +
  ylab("Landscape diversity") +
  xlab("Time (days)") +
  ylab("Landscape diversity") +
  labs(colour = "Corridor length", fill = "Corridor length") +
  theme_classic()+
  theme(legend.position = "top",
        aspect.ratio = 1,
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7)) +
  scale_colour_manual(values = c("#3B528BFF", "#C7E020FF")) 

### Combine temperature and corridors time plots and export ###
time_plots <- ggpubr::ggarrange(temperature_time, corridor_time, labels = "auto", label.x = 0.05, label.y = 0.72)

ggsave("Results/Figure2.tiff", time_plots, units = "in", width = 10, height = 10)
