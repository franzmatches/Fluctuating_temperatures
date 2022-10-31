### Figure 2 ###

## Preamble
library(ggh4x)
library(tidyverse) # data wrangling and plotting

source("Code/data_preparation.R")
palette <- c("#440154FF", "#3B528BFF","#21908CFF", "#5DC863FF")

### Species diversity across temperature regimes model ###
modTemp <- glmmTMB::glmmTMB(Sh_div~ splines::ns(NumDays,3) + Temp_Regime + 
                            splines::ns(NumDays,3):Temp_Regime +
                            ar1(as.factor(NumDays) + 0 | Replicate) + (1|Replicate), 
                          data = data_pooled,family = "gaussian",REML=F)

gg.df_TEMP <- data.frame(preds = predict(modTemp, newdata = data_pooled,type = "response",re.form = NA),
                         upr = predict(modTemp, newdata = data_pooled ,se.fit = T,type = "response", re.form = NA)$fit +
                           predict(modTemp, newdata = data_pooled,se.fit = T,type = "response", re.form = NA)$se.fit ,
                         lwr = predict(modTemp, newdata = data_pooled ,se.fit = T ,type = "response", re.form = NA)$fit -
                           predict(modTemp, newdata = data_pooled ,se.fit = T ,type = "response", re.form = NA)$se.fit , 
                         data_pooled)

temperature_time <- ggplot(gg.df_TEMP, aes(x =NumDays, y= Sh_div,group = Temp_Regime,
                                           col =Temp_Regime))+
  geom_point(alpha = 1, pch = 1)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour = NA, alpha=0.1)+
  xlab("Time (Days)") +
  ylab("Shannon diversity") +
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
mod_corr<-glmmTMB::glmmTMB(Sh_div~ splines::ns(NumDays,3) + Corridor + 
                             splines::ns(NumDays,3):Corridor +
                             ar1(as.factor(NumDays) + 0 | Replicate) + (1|Replicate), 
                           data = data_pooled,family = "gaussian",REML=F)


gg.df_Cor <- data.frame(preds = predict(mod_corr, newdata = data_pooled,type = "response", re.form = NA),
                        upr = predict(mod_corr, newdata = data_pooled ,se.fit = T,type = "response", re.form = NA)$fit +
                          predict(mod_corr, newdata = data_pooled,se.fit = T,type = "response", re.form = NA)$se.fit ,
                        lwr = predict(mod_corr, newdata = data_pooled ,se.fit = T ,type = "response", re.form = NA)$fit -
                          predict(mod_corr, newdata = data_pooled ,se.fit = T ,type = "response", re.form = NA)$se.fit , 
                        data_pooled)

corridor_time <- ggplot(gg.df_Cor, aes(x =NumDays, y= Sh_div,group = Corridor,
                                       colour =Corridor))+
  geom_point(alpha = 1, pch = 1)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr,
                  group = Corridor),colour = NA, alpha=0.1)+
  xlab("Time (Days)") +
  ylab("Shannon diversity") +
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
