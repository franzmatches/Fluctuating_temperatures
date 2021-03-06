### Figure 4 ###

## Preamble
library(ggh4x)
library(tidyverse) # data wrangling and plotting

source("Code/data_preparation.R")
palette <- c("#440154FF", "#3B528BFF","#21908CFF", "#5DC863FF")

colnames(data_pooled)[6:12] <- c("Blepharisma","Didinium","Colpidium","Homolazoon","Bursaria","Spirostomum","Paramecium") #rename species for downstream plotting
data_pooled_lastday <- data_pooled %>%
  filter(NumDays ==max(NumDays)) #filter to last day of experiment

data_pooled_half <- data_pooled %>%
  filter(NumDays == max(NumDays)/2) #filter to halfway through experiment

### Community composition across treatments at end of experiment ###
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

### Community composition across treatments at end of experiment ###
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
  geom_point(aes(alpha = Corridor,fill=Temp_Regime,col=Temp_Regime),shape = 21,size = 1.5) + # points for Long Corridors
  geom_point(aes(col=Temp_Regime),shape = 21,size = 1.5) + # points for Short Corridors
  geom_point(aes(x = PC1.centroid, y= PC2.centroid,fill=Temp_Regime,col=Temp_Regime),shape = 8,size = 1.5) + # points for Temp_Regime centroids
  ggpubr::stat_chull(aes(fill=Temp_Regime,col=Temp_Regime), alpha = 0.1, 
                     geom = "polygon")+ # Temp_Regime convex hulls
  xlab("PC1 (26.6% var explained)") + ylab("PC2 (22.1% var explained)")+
  scale_alpha_manual(values = c(1,0),
                     guide = guide_legend(override.aes = list(alpha = 1,fill = c("black",NA))))+
  geom_segment(data = pca_half_vars,aes(x = 0, y= 0, xend = PC1,yend=PC2),
               arrow = arrow(length=unit(0.15,"cm")),alpha=0.7)+ # add species contribution arrows
  geom_text(data= pca_half_vars, aes(x=PC1, y=PC2, label=species), 
            size = 3,color="black",
            hjust = c(-0.1, 0.8, -0.4, -0.2, 1, 0.8, 1), #specify positions for each label separately because they are still overlapping with ggrepel
            vjust = c(1.5, 1.5, 1.3, 0.6, -0.3, 2, -0.4))+
  xlim(-5.5,5.5)+ylim(-7.0,7.0)+ #square plot area for symmetry and unbiased interpretation
  theme_classic()+
  theme(legend.position = c(0.22, 0.75),
        legend.key.size = unit(0.4, 'cm'),
        legend.key.width = unit(0.4, 'cm'),
        legend.key.height = unit(0.4, 'cm'),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        aspect.ratio = 1) +
  labs(colour = "Temperature regime", fill = "Temperature regime", alpha = "Corridor length")+
  scale_colour_manual(values = palette, labels = c("Constant", "Fluctuating asynchronous", "Fluctuating synchronous", "Static difference"))+
  scale_fill_manual(values = palette, labels = c("Constant", "Fluctuating asynchronous", "Fluctuating synchronous", "Static difference"))

## Combine plots ##
PCA_plots <- ggpubr::ggarrange(pca_half.plot, pca_lastday.plot, ncol = 1, align = "hv", labels = "auto", label.x = 0.25, label.y = 0.95)
ggsave("Results/Figure4.tiff", PCA_plots, units = "in", width = 10, height = 10)
