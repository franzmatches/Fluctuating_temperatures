#########################################################################################
### Community ###
#########################################################################################

### Preamble ###
library(tidyverse) # data wrangling and plotting
library(vegan) # community ecology metrics and analysis
library(ggrepel) # non-overlapping labels
library(ggpubr) # multi-panelled plotting
library(FactoMineR)
library(factoextra)

source("Code/data_preparation.R") # load raw data

colnames(data_pooled)[6:12] <- c("Blepharisma","Didinium","Colpidium","Homolazoon","Bursaria","Spirostomum","Paramecium") #rename species for downstream plotting
data_pooled_lastday <- data_pooled %>%
  filter(NumDays ==max(NumDays)) #filter to last day of experiment

data_pooled_half <- data_pooled %>%
  filter(NumDays == max(NumDays)/2) #filter to halfway through experiment

#----------------------------------------------------------------------------------------
## VISUALISE COMMUNITY AT THE END OF THE EXPERIMENT ## 
#----------------------------------------------------------------------------------------

week_four_microcosm_PCA <- prcomp(data_pooled_lastday[,c(6,8:12)],scale. = T) #estimate PCs using base R prcomp (scaled and centered)
#Didinium dropped due to all zeroes across treatments and inability to scale when variance is constant

summary(week_four_microcosm_PCA) # PCs 1 & 2  contribute ~34 and 19% variation respectively

pca_lastday_df <- data.frame(PC1 = week_four_microcosm_PCA$x[,1], #extract replicate positions and first two PC axes
                             PC2 = week_four_microcosm_PCA$x[,2], 
                             data_pooled_lastday) %>% # combine with original dataset
  left_join(aggregate(cbind(PC1.centroid = week_four_microcosm_PCA$x[,1],PC2.centroid = week_four_microcosm_PCA$x[,2])~Temp_Regime, 
                      data =data_pooled_lastday, mean)) #extract centroid positions (mean) for plotting (if desired)


pca_lastday_vars <- data.frame(PC1 = week_four_microcosm_PCA$rotation[,1]*5, #extract species contributions and multiply by 5 for prominent arrows
                               PC2 =  week_four_microcosm_PCA$rotation[,2]*5,
                               species = rownames(week_four_microcosm_PCA$rotation)) # label contributions with associated species

palette <- c("#440154FF", "#3B528BFF","#21908CFF", "#5DC863FF")

#recreate fviz_pca_biplot as raw ggplot. Options provided for either 95% multinormal ellipses or convex hull
pca_lastday.plot <- ggplot(pca_lastday_df, aes(x = PC1,y=PC2))  + 
  geom_vline(xintercept = 0,linetype = "dashed",col="black",alpha=0.8)+ # lines denoting community space quadrants
  geom_hline(yintercept = 0,linetype = "dashed",col="black",alpha=0.8)+
  geom_point(aes(alpha = Corridor,fill=Temp_Regime,col=Temp_Regime),shape = 21,size = 1.5) + # points for Long Corridors
  geom_point(aes(col=Temp_Regime),shape = 21,size = 1.5) + # points for Short Corridors
  geom_point(aes(x = PC1.centroid, y= PC2.centroid,fill=Temp_Regime,col=Temp_Regime),shape = 8,size = 1.5) + # points for Temp_Regime centroids
  ggpubr::stat_chull(aes(fill=Temp_Regime,col=Temp_Regime), alpha = 0.1, 
                     geom = "polygon")+ # Temp_Regime convex hulls
  xlab("PC1 (34.7% var explained)") + ylab("PC2 (19.1% var explained)")+
  scale_alpha_manual(values = c(1,0),
                     guide = guide_legend(override.aes = list(alpha = 1,fill = c("black",NA))))+ # legend customisation of Long vs Short corridors
  geom_segment(data = pca_lastday_vars,aes(x = 0, y= 0, xend = PC1,yend=PC2),
               arrow = arrow(length=unit(0.15,"cm")),alpha=0.7)+ # add species contribution arrows
  geom_text(data= pca_lastday_vars, aes(x=PC1, y=PC2, label=species), 
            size = 3,color="black",
            hjust = c(-0.1, -0.5, 0.2, 1.3, 1.2, 1), #specify positions for each label separately because they are still overlapping with ggrepel
            vjust = c(1, 0.5, 2, 1, -1.2, 2))+
  xlim(-6.5,6.5)+ylim(-5,5)+ #symmetrical plot area for unbiased interpretation
  theme_classic() +
  theme(legend.position = "none",
        aspect.ratio = 1) +
  labs(colour = "Temperature regime", fill = "Temperature regime", alpha = "Corridor length")+
  scale_colour_manual(values = palette, labels = c("Constant", "Fluctuating asynchronous", "Fluctuating synchronous", "Static difference"))+
  scale_fill_manual(values = palette, labels = c("Constant", "Fluctuating asynchronous", "Fluctuating synchronous", "Static difference"))

## ASSESS EFFECT OF TREATMENT ON FINAL COMMUNITY COMPOSITION ##
pooled_lastday_dist <- vegan::vegdist(data_pooled_lastday[,c(6:12)],method = "bray") #estimate bray-curtis dissimilarity due high proportion of zeroes

#perform PERMANOVA and PERMDIST.
vegan::adonis2(pooled_lastday_dist~data_pooled_lastday$Temp_Regime + data_pooled_lastday$Corridor,
               data = data_pooled_lastday,by="margin") # no interaction included per other analyses
#adonis2 preferred to adonis due to improved customisation and clarity. 
#by = "margin" is marginal effect of each predictor rather than additional effect that
#predictor has over previous. Therefore order of predictors supplied to model irrelevant

#PERMANOVA pairwise comparisons
library(pairwiseAdonis)
#last day
pairwise.adonis(data_pooled_lastday[,c(6,8:12)], data_pooled_lastday$Temp_Regime)
#half way
pairwise.adonis(data_pooled_half[,c(6:12)], data_pooled_half$Temp_Regime)

#!!!DEFINING `TYPE` AS CENTROID/MEDIAN ALTERS INTERPRETATION!!!# Needs discussing - decided to stick with centroid
vegan::permutest(vegan::betadisper(pooled_lastday_dist,data_pooled_lastday$Temp_Regime,type = "centroid"),pairwise = T) #compare variance of Temp_Regime groups 
#PERMDIST with pairwise comparisons. permutest = vegan specific function with greater customisation
#dispersions are heterogeneous due Constant:Fluctuating_Synchro
#TukeyHSD(vegan::betadisper(pooled_lastday_dist,data_pooled_lastday$Temp_Regime,type = "centroid")) #agrees

vegan::permutest(vegan::betadisper(pooled_lastday_dist,data_pooled_lastday$Corridor,type = "centroid")) #compare variance of Corridor groups
#dispersions are homogeneous

#----------------------------------------------------------------------------------------
## VISUALISE COMMUNITY HALFWAY THROUGH THE EXPERIMENT ## 
#----------------------------------------------------------------------------------------

#week_two_microcosm_PCA <- FactoMineR::PCA(data_pooled_half[,c(6:12)],scale.unit = TRUE, graph = FALSE)
week_two_microcosm_PCA <- prcomp(data_pooled_half[,c(6:12)],scale. = T)

summary(week_two_microcosm_PCA) # PCs 1 & 2  contribute ~27 and 22% variation respectively

pca_half_df <- data.frame(PC1 = week_two_microcosm_PCA$x[,1],
                          PC2 = week_two_microcosm_PCA$x[,2], 
                          data_pooled_half) %>% # combine with original dataset
  left_join(aggregate(cbind(PC1.centroid = week_two_microcosm_PCA$x[,1],PC2.centroid = week_two_microcosm_PCA$x[,2])~Temp_Regime, 
                      data =data_pooled_half, mean)) #extract centroid positions (mean) for plotting (if desired)

pca_half_vars <- data.frame(PC1 = week_two_microcosm_PCA$rotation[,1]*5, #extract species contributions and multiply by 5 for prominent arrows
                            PC2 =  week_two_microcosm_PCA$rotation[,2]*5,
                            species = rownames(week_two_microcosm_PCA$rotation))

pca_half.plot <- ggplot(pca_half_df, aes(x = PC1,y=PC2))  + 
  geom_vline(xintercept = 0,linetype = "dashed",col="black",alpha=0.8)+
  geom_hline(yintercept = 0,linetype = "dashed",col="black",alpha=0.8)+
  #geom_point(aes(alpha = Corridor,fill=Temp_Regime,col=Temp_Regime),shape = 21,size = 3, position=position_jitter(0.8,seed=123)) + # points for Long Corridors
  #geom_point(aes(col=Temp_Regime),shape = 21,size = 3, position=position_jitter(0.8,seed=123)) + # points for Short Corridors
  geom_point(aes(alpha = Corridor,fill=Temp_Regime,col=Temp_Regime),shape = 21,size = 3) + # points for Long Corridors
  geom_point(aes(col=Temp_Regime),shape = 21,size = 3) + # points for Short Corridors
  stat_ellipse(aes(fill=Temp_Regime,col=Temp_Regime), alpha=.2,type='t',size =0.3, geom="polygon",level = 0.95)+ # 95% multinormal ellipses 
  #geom_point(aes(x = PC1.centroid, y= PC2.centroid,fill=Temp_Regime,col=Temp_Regime),shape = 8,size = 5) + # points for Temp_Regime centroids
  #ggpubr::stat_chull(aes(fill=Temp_Regime,col=Temp_Regime), alpha = 0.1, 
  #                   geom = "polygon")+ # Temp_Regime convex hulls
  xlab("PC1 (26.6% var explained)") + ylab("PC2 (22.1% var explained)")+
  scale_alpha_manual(values = c(1,0),
                     guide = guide_legend(override.aes = list(alpha = 1,fill = c("black",NA))))+
  geom_segment(data = pca_half_vars,aes(x = 0, y= 0, xend = PC1,yend=PC2),
               arrow = arrow(length=unit(0.15,"cm")),alpha=0.7)+ # add species contribution arrows
  #geom_text(data= pca_half_vars, aes(x=PC1, y=PC2, label=species), size = 3, vjust =1,color="black")+ # label arrows
  ggrepel::geom_text_repel(data= pca_half_vars, aes(x=PC1, y=PC2, label=species), size = 3,color="black",direction = "y")+ # label arrows
  xlim(-5.5,5.5)+ylim(-7.0,7.0)+ #square plot area for symmetry and unbiased interpretation
  theme_bw()

## ASSESS EFFECT OF TREATMENT ON INTERMEDIATE COMMUNITY COMPOSITION ##
pooled_half_dist <- vegan::vegdist(data_pooled_half[,c(6:12)],method = "bray")
#perform PERMANOVA and PERMDIST
vegan::adonis2(pooled_half_dist~data_pooled_half$Temp_Regime + data_pooled_half$Corridor,
               data = data_pooled_half, method = "bray",by="margin") # no interaction per other analyses
vegan::permutest(vegan::betadisper(pooled_half_dist,data_pooled_half$Temp_Regime,type = "centroid"),pairwise = T) #PERMDIST with pairwise comparisons 
#TukeyHSD(vegan::betadisper(pooled_half_dist,data_pooled_half$Temp_Regime,type = "centroid")) #agrees

vegan::permutest(vegan::betadisper(pooled_half_dist,data_pooled_half$Corridor,type = "centroid"),pairwise = T) 
#TukeyHSD(vegan::betadisper(pooled_half_dist,data_pooled_half$Corridor,type = "centroid"))
#dispersions are hetergeneous


#----------------------------------------------------------------------------------------
## Combined 14 and 28 day PCA ##
ggpubr::ggarrange(pca_half.plot,pca_lastday.plot,nrow = 2,ncol = 1,common.legend = T, labels = c("a","b"))
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
## VISUALISE AND ASSESS CHANGE OVER TIME ## 
#Principal response curves (DOI 10.1007/s10661-008-0314-6)
#----------------------------------------------------------------------------------------

#generate prcs for corridor and temperature regime separately, and then for an interaction factor 
#(prcs can currently only be performed on a single treatment). prc compares influence of each treatment
#relative to a control over time by conditioning the fit on the time variable (i.e. NumDays)
corridor.prc <- vegan::prc(log1p(data_pooled[,c(6:12)]),treatment = data_pooled$Corridor, as.factor(data_pooled$NumDays)) #control = Long corridors
temp_regime.prc <- vegan::prc(log1p(data_pooled[,c(6:12)]),treatment = data_pooled$Temp_Regime, as.factor(data_pooled$NumDays)) #control = Constant temp_regime
interaction.prc <-  vegan::prc(log1p(data_pooled[,c(6:12)]),treatment = interaction(data_pooled$Corridor,data_pooled$Temp_Regime), as.factor(data_pooled$NumDays)) #control = Long.Constant
#multiple treatment factors aren't possible in standard prc, therefore explore interactions

summary_temp_regime1 <- summary(temp_regime.prc, scaling = "symmetric", axis = 1,
                                correlation = F) #extract first component 
summary_temp_regime2 <- summary(temp_regime.prc, scaling = "symmetric", axis = 2,
                                correlation = F) #extract second component 
summary_corridor1 <- summary(corridor.prc, scaling = "symmetric", axis = 1,
                             correlation = F)
summary_corridor2 <- summary(corridor.prc, scaling = "symmetric", axis = 2,
                             correlation = F)
summary_interaction1 <- summary(interaction.prc, scaling = "symmetric", axis = 1,
                                correlation = F)
summary_interaction2 <- summary(interaction.prc, scaling = "symmetric", axis = 2,
                                correlation = F)

anova(corridor.prc) # differences identified using an ANOVA but no post hoc tests available
anova(temp_regime.prc) # differences identified
anova(interaction.prc) # differences identified

# visualise PRCs
df_pc1 <- data.frame(NumDays = unique(data_pooled$NumDays), Short = t(coef(summary_corridor1)), Long = 0, Constant = 0) %>% 
  cbind(t(coef(summary_temp_regime1))) %>% #merge time data, control effects and estimated effects from prc
  pivot_longer(-NumDays,names_to = "treatment",values_to = "effect") %>% # pivot to long format
  mutate(Treatment_category = ifelse(treatment %in% c("Short","Long"),"Corridor","Temp_Regime"))%>% # classify into treatments
  mutate(PC = "PC1") # classify associated principal component

df_pc2 <- data.frame(NumDays = unique(data_pooled$NumDays), Short = t(coef(summary_corridor2)), Long = 0, Constant = 0) %>% 
  cbind(t(coef(summary_temp_regime2))) %>%
  pivot_longer(-NumDays,names_to = "treatment",values_to = "effect") %>%
  mutate(Treatment_category = ifelse(treatment %in% c("Short","Long"),"Corridor","Temp_Regime"))%>%
  mutate(PC = "PC2")

df_fin <- rbind(df_pc1,df_pc2)

#visualise estimated effect differences of treatments vs control treatment (Long corridor for corridor and Constant for temperature)
ggplot(df_fin,aes(x=NumDays,y=effect)) + 
  geom_point(aes(col=treatment))+ # points and...
  geom_line(aes(col=treatment)) + # lines of effects through time
  facet_grid(PC~Treatment_category) + #combine both treatment types in to single figure
  scale_color_manual(values = scales::hue_pal()(6),name = "Treatment") +
  theme_bw() + xlab("Day of experiment") + ylab("Effect") +
  geom_rug(data = rbind(data.frame(spp = names(summary_corridor1$sp),effect = summary_corridor1$sp*0.1, PC = "PC1",Treatment_category = "Corridor"),
                        data.frame(spp = names(summary_corridor1$sp),effect = summary_corridor1$sp*0.1, PC = "PC2",Treatment_category = "Corridor"),
                        data.frame(spp = names(summary_temp_regime1$sp),effect = summary_temp_regime1$sp*0.1, PC = "PC1",Treatment_category = "Temp_Regime"),
                        data.frame(spp = names(summary_temp_regime2$sp),effect = summary_temp_regime2$sp*0.1, PC = "PC2",Treatment_category = "Temp_Regime")),
           sides = "r",outside = F, #geom_rug plots influence strength of individual species on second y axis
           mapping = aes_string(group = NULL, x = NULL,
                                colour = NULL, linetype = NULL)) +
  scale_y_continuous("Effect", sec.axis = sec_axis(trans = ~ .*10, name = "Species influence")) + # plot and standardise second y axis as magnitude of effects too large to fit on original (left hand) scale
  ggrepel::geom_text_repel(data= rbind(data.frame(spp = names(summary_corridor1$sp),effect = summary_corridor1$sp*0.1, PC = "PC1",Treatment_category = "Corridor"),
                                       data.frame(spp = names(summary_corridor1$sp),effect = summary_corridor1$sp*0.1, PC = "PC2",Treatment_category = "Corridor"),
                                       data.frame(spp = names(summary_temp_regime1$sp),effect = summary_temp_regime1$sp*0.1, PC = "PC1",Treatment_category = "Temp_Regime"),
                                       data.frame(spp = names(summary_temp_regime2$sp),effect = summary_temp_regime2$sp*0.1, PC = "PC2",Treatment_category = "Temp_Regime")),
                           mapping=aes(x=32,y=effect,label=spp),hjust=-0.5, size=3, direction = "y",box.padding = 0.3,min.segment.length = 0.1,vjust=0.55) + # label species and ensure minimal overlap
  xlim(0,32) # extend y axis to fit labels

#visualise estimated effect differences of interacting treatments vs control treatment (Long corridor.Constant temperature)
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
                           mapping=aes(x=30,y=effect,label=spp),hjust=-1, size=3, direction = "y",box.padding = 0.1,min.segment.length = 0.05,vjust = 0.55) 

#end