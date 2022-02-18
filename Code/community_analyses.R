#########################################################################################
### Community ###
#########################################################################################

### Preamble ###
library(tidyverse) # data wrangling and plotting
library(FactoMineR)
library(factoextra)

source("Code/data_preparation.R")

colnames(data_pooled)[6:12] <- c("Blepharisma","Didinium","Colpidium","Homolazoon","Bursaria","Spirostomum","Paramecium")
data_pooled_lastday <- data_pooled %>%
  filter(NumDays ==max(NumDays)) 

data_pooled_half <- data_pooled %>%
  filter(NumDays == max(NumDays)/2) 

#----------------------------------------------------------------------------------------
## VISUALISE COMMUNITY AT THE END OF THE EXPERIMENT ## 
#----------------------------------------------------------------------------------------

#week_four_microcosm_PCA <- FactoMineR::PCA(data_pooled_lastday[,c(6:12)],scale.unit = TRUE, graph = FALSE)
week_four_microcosm_PCA <- prcomp(data_pooled_lastday[,c(6,8:12)],scale. = T)

summary(week_four_microcosm_PCA) # PCs 1 & 2  contribute ~34 and 19% variation respectively

pca_lastday_df <- data.frame(PC1 = week_four_microcosm_PCA$x[,1],
                             PC2 = week_four_microcosm_PCA$x[,2], 
                             data_pooled_lastday)

pca_lastday.plot <- ggplot(pca_lastday_df, aes(x = PC1,y=PC2,col=Temp_Regime))  + 
  geom_point(aes(shape = Corridor),size = 3, position=position_jitter(0.8)) +
  stat_ellipse(aes(fill=Temp_Regime), alpha=.2,type='t',size =0.3, geom="polygon")+
  xlab("PC 1 (34.7% var explained)") + ylab("PC2 (19.1% var explained)")+
  theme_bw()

# fviz_pca_biplot(week_four_microcosm_PCA,
#                 col.var = "black",
#                 label = "var",
#                 alpha.ind = 0.1,#make these points invisible to have diff shapes for points
#                 labelsize = 2,
#                 repel = TRUE,
#                 geom_ind = "points",
#                 col.ind = data_pooled_lastday$Temp_Regime,
#                 palette = c("#00AFBB", "#E7B800", "#FC4E07", "#0073C2FF"),
#                 mean.point = FALSE,
#                 addEllipses = TRUE,
#                 title = NULL,
#                 ellipse.type = "confidence",
#                 legend.title = "Temperature treatment",
#                 geom.ind = "point")+
#   #add different shape points for long vs. short corridors (can't do this with fviz)
#   geom_point(aes(shape = data_pooled_lastday$Corridor,
#                  colour = data_pooled_lastday$Temp_Regime)) +
#   scale_shape_manual(values = c(3,4,15,17,19,16))+
#   theme(legend.position = "none",
#         aspect.ratio = 1)+
#   labs(x = "PC1", y = "PC2")

adonis(data_pooled_lastday[,c(6:12)]~data_pooled_lastday$Temp_Regime + data_pooled_lastday$Corridor,
       data = data_pooled_lastday, method = "bray") # no interaction per other analyses

#----------------------------------------------------------------------------------------
## VISUALISE COMMUNITY HALFWAY THROUGH THE EXPERIMENT ## 
#----------------------------------------------------------------------------------------

#week_four_microcosm_PCA <- FactoMineR::PCA(data_pooled_lastday[,c(6:12)],scale.unit = TRUE, graph = FALSE)
week_two_microcosm_PCA <- prcomp(data_pooled_half[,c(6:12)],scale. = T)

summary(week_two_microcosm_PCA) # PCs 1 & 2  contribute ~27 and 22% variation respectively

pca_half_df <- data.frame(PC1 = week_two_microcosm_PCA$x[,1],
                             PC2 = week_two_microcosm_PCA$x[,2], 
                             data_pooled_half)

pca_half.plot <- ggplot(pca_half_df, aes(x = PC1,y=PC2,col=Temp_Regime))  + 
  geom_point(aes(shape = Corridor),size = 3, position=position_jitter(0.8)) +
  stat_ellipse(aes(fill=Temp_Regime), alpha=.2,type='t',size =0.3, geom="polygon")+
  xlab("PC 1 (26.6% var explained)") + ylab("PC2 (22.1% var explained)")+
  theme_bw()

# fviz_pca_biplot(week_two_microcosm_PCA,
#                 col.var = "black",
#                 label = "var",
#                 alpha.ind = 0.1,#make these points invisible to have diff shapes for points
#                 labelsize = 2,
#                 repel = TRUE,
#                 geom_ind = "points",
#                 col.ind = data_pooled_half$Temp_Regime,
#                 palette = c("#00AFBB", "#E7B800", "#FC4E07", "#0073C2FF"),
#                 mean.point = FALSE,
#                 addEllipses = TRUE,
#                 title = NULL,
#                 ellipse.type = "confidence",
#                 legend.title = "Temperature treatment",
#                 geom.ind = "point")+
#   #add different shape points for long vs. short corridors (can't do this with fviz)
#   geom_point(aes(shape = data_pooled_half$Corridor,
#                  colour = data_pooled_half$Temp_Regime)) +
#   scale_shape_manual(values = c(3,4,15,17,19,16))+
#   theme(legend.position = "none",
#         aspect.ratio = 1)+
#   labs(x = "PC1", y = "PC2")

adonis(data_pooled_half[,c(6:12)]~data_pooled_half$Temp_Regime + data_pooled_half$Corridor,
       data = data_pooled_half, method = "bray") # no interaction per other analyses

#----------------------------------------------------------------------------------------
## Combined 14 and 28 day PCA ##
ggpubr::ggarrange(pca_half.plot,pca_lastday.plot,nrow = 2,ncol = 1,common.legend = T, labels = c("a","b"))
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
## VISUALISE AND ASSESS CHANGE OVER TIME ## 
#Principal response curves (DOI 10.1007/s10661-008-0314-6)
#----------------------------------------------------------------------------------------

#generate prcs for corridor and temperature regime separately, and then for an interaction factor 
#(prcs can currently only be performed on a single treatment)
corridor.prc <- vegan::prc(log1p(data_pooled[,c(6:12)]),treatment = data_pooled$Corridor, as.factor(data_pooled$NumDays))
temp_regime.prc <- vegan::prc(log1p(data_pooled[,c(6:12)]),treatment = data_pooled$Temp_Regime, as.factor(data_pooled$NumDays))
interaction.prc <-  vegan::prc(log1p(data_pooled[,c(6:12)]),treatment = interaction(data_pooled$Corridor,data_pooled$Temp_Regime), as.factor(data_pooled$NumDays))

summary_temp_regime1 <- summary(temp_regime.prc, scaling = "symmetric", axis = 1,
                              correlation = F) #extract first component 
summary_temp_regime2 <- summary(temp_regime.prc, scaling = "symmetric", axis = 2,
                               correlation = F) #extract first component 
summary_corridor1 <- summary(corridor.prc, scaling = "symmetric", axis = 1,
                              correlation = F)
summary_corridor2 <- summary(corridor.prc, scaling = "symmetric", axis = 2,
                             correlation = F)
summary_interaction1 <- summary(interaction.prc, scaling = "symmetric", axis = 1,
                             correlation = F)
summary_interaction2 <- summary(interaction.prc, scaling = "symmetric", axis = 2,
                             correlation = F)

anova(corridor.prc) # differences identified using an ANOVA but no post hoc tests available
permutest(corridor.prc)
anova(temp_regime.prc) # differences identified
anova(interaction.prc) # differences identified

# visualise PRCs
df_pc1 <- data.frame(NumDays = unique(data_pooled$NumDays), Short = t(coef(summary_corridor1)), Long = 0, Constant = 0) %>% 
  cbind(t(coef(summary_temp_regime1))) %>%
  pivot_longer(-NumDays,names_to = "treatment",values_to = "effect") %>%
  mutate(Treatment_category = ifelse(treatment %in% c("Short","Long"),"Corridor","Temp_Regime"))%>%
  mutate(PC = "PC1")
  
df_pc2 <- data.frame(NumDays = unique(data_pooled$NumDays), Short = t(coef(summary_corridor2)), Long = 0, Constant = 0) %>% 
  cbind(t(coef(summary_temp_regime2))) %>%
  pivot_longer(-NumDays,names_to = "treatment",values_to = "effect") %>%
  mutate(Treatment_category = ifelse(treatment %in% c("Short","Long"),"Corridor","Temp_Regime"))%>%
  mutate(PC = "PC2")

df_fin <- rbind(df_pc1,df_pc2)

#visualise differences of treatments vs control treatment (Long corridor for corridor and Constant for temperature)
ggplot(df_fin,aes(x=NumDays,y=effect)) + 
  #geom_line(aes(linetype = Treatment_category)) + facet_wrap(~PC) +
  geom_point(aes(col=treatment))+
  geom_line(aes(col=treatment)) + facet_grid(PC~Treatment_category) + #combine both treatment types in to single figure
  scale_color_manual(values = scales::hue_pal()(6),name = "Treatment") +
  theme_bw() + xlab("Day of experiment") + ylab("Effect") +
  geom_rug(data = rbind(data.frame(spp = names(summary_corridor1$sp),effect = summary_corridor1$sp*0.1, PC = "PC1",Treatment_category = "Corridor"),
                       data.frame(spp = names(summary_corridor1$sp),effect = summary_corridor1$sp*0.1, PC = "PC2",Treatment_category = "Corridor"),
                       data.frame(spp = names(summary_temp_regime1$sp),effect = summary_temp_regime1$sp*0.1, PC = "PC1",Treatment_category = "Temp_Regime"),
                       data.frame(spp = names(summary_temp_regime2$sp),effect = summary_temp_regime2$sp*0.1, PC = "PC2",Treatment_category = "Temp_Regime")),
          sides = "r",outside = F, #geom_rug plots influence strength of species on second y axis
          mapping = aes_string(group = NULL, x = NULL,
                               colour = NULL, linetype = NULL)) +
  scale_y_continuous("Effect", sec.axis = sec_axis(trans = ~ .*10, name = "Species influence")) + #standardise second y axis
  ggrepel::geom_text_repel(data= rbind(data.frame(spp = names(summary_corridor1$sp),effect = summary_corridor1$sp*0.1, PC = "PC1",Treatment_category = "Corridor"),
                                       data.frame(spp = names(summary_corridor1$sp),effect = summary_corridor1$sp*0.1, PC = "PC2",Treatment_category = "Corridor"),
                                       data.frame(spp = names(summary_temp_regime1$sp),effect = summary_temp_regime1$sp*0.1, PC = "PC1",Treatment_category = "Temp_Regime"),
                                       data.frame(spp = names(summary_temp_regime2$sp),effect = summary_temp_regime2$sp*0.1, PC = "PC2",Treatment_category = "Temp_Regime")),
                           mapping=aes(x=30,y=effect,label=spp),hjust=0, size=3, direction = "y",box.padding = 0.3,min.segment.length = 0.1,vjust=0.55) #label species

df_interaction_pc1 <- data.frame(NumDays = unique(data_pooled$NumDays), Long.Constant =0) %>% 
  cbind(t(coef(summary_interaction1))) %>%
  pivot_longer(-NumDays,names_to = "treatment",values_to = "effect") %>%
  mutate(PC = "PC1")

df_interaction_pc2 <- data.frame(NumDays = unique(data_pooled$NumDays), Long.Constant =0) %>% 
  cbind(t(coef(summary_interaction2))) %>%
  pivot_longer(-NumDays,names_to = "treatment",values_to = "effect") %>%
  mutate(PC = "PC2")

df_interaction_fin <- rbind(df_interaction_pc1,df_interaction_pc2)

ggplot(df_interaction_fin,aes(x=NumDays,y=effect)) + 
  #geom_line(aes(linetype = Treatment_category)) + facet_wrap(~PC) +
  geom_point(aes(col=treatment))+
  geom_line(aes(col=treatment)) + facet_wrap(~PC,strip.position = "right",dir="v") + 
  scale_color_manual(values = scales::hue_pal()(8),name = "Treatment") +
  theme_bw() + xlab("Day of experiment") + ylab("Effect") +
  geom_rug(data = rbind(data.frame(spp = names(summary_interaction1$sp),effect = summary_interaction1$sp*0.05, PC = "PC1"),
                        data.frame(spp = names(summary_interaction2$sp),effect = summary_interaction2$sp*0.05, PC = "PC2")),
           sides = "r",outside = F,
           mapping = aes_string(group = NULL, x = NULL,
                                colour = NULL, linetype = NULL)) +
  scale_y_continuous("Effect", sec.axis = sec_axis(trans = ~ .*20, name = "Species influence")) + 
  ggrepel::geom_text_repel(data= rbind(data.frame(spp = names(summary_interaction1$sp),effect = summary_interaction1$sp*0.05, PC = "PC1"),
                        data.frame(spp = names(summary_interaction2$sp),effect = summary_interaction2$sp*0.05, PC = "PC2")),
            mapping=aes(x=30,y=effect,label=spp),hjust=0, size=3, direction = "y",box.padding = 0.2,min.segment.length = 0.1,vjust = 0.1) 
