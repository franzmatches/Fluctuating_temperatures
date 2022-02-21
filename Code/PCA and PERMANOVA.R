#load packages
library(tidyverse)
library(FactoMineR)
library(factoextra)

#patch-level data through time
patch_dataset <- read.csv("Data/fluctuation_time_abundances.csv",
                          header = T,
                          row.names = 1)

#group together to microcosm level
microcosm_dataset <- patch_dataset %>%
  #species presence for each corridor length, temp treatment, replicate, day
  group_by(corridor_length, temp_treatment, replicate, day) %>%
  #want the sum of column for each individual species
  summarise(blepharisma = sum(blepharisma),
            colpidium = sum(colpidium),
            didinium = sum(didinium),
            homalozoon = sum(homalozoon),
            paramecium = sum(paramecium),
            spirostomum = sum(spirostomum),
            bursaria = sum(bursaria))

### Subset by time - final day ####
#subset final day
week_four_microcosm_dataset <- 
  microcosm_dataset %>% filter(day == 28)
#remove didinium because they were completely absent
week_four_microcosm_dataset$didinium <- NULL

#run the PCA
week_four_microcosm_PCA <- PCA(week_four_microcosm_dataset[,c(5:10)], graph = FALSE)

#examine eigenvalues to determine number of PCs
week_four_eig_vals <- get_eigenvalue(week_four_microcosm_PCA)
week_four_eig_vals #PCs 1, 2 and 3 have eigenvalues greater than 1 and account for ~72% variance

#extract the results
week_four_vars <- get_pca_var(week_four_microcosm_PCA)

#correlation between variables and PCs
week_four_correlations <- week_four_vars$coord

#quality of representation
week_four_representation <- week_four_vars$cos2

#contribution of variables to PCs
week_four_contributions <- week_four_vars$contrib

#biplot of microcosms and species with confidence ellipses

week_four_biplot <- fviz_pca_biplot(week_four_microcosm_PCA,
                                    col.var = "black",
                                    label = "var",
                                    alpha.ind = 0.1,#make these points invisible to have diff shapes for points
                                    labelsize = 2,
                                    repel = TRUE,
                                    geom_ind = "points",
                                    col.ind = week_four_microcosm_dataset$temp_treatment,
                                    palette = c("#00AFBB", "#E7B800", "#FC4E07", "#0073C2FF"),
                                    mean.point = FALSE,
                                    addEllipses = TRUE,
                                    title = NULL,
                                    ellipse.type = "confidence",
                                    legend.title = "Temperature treatment",
                                    geom.ind = "point")+
  #add different shape points for long vs. short corridors (can't do this with fviz)
  geom_point(aes(shape = week_four_microcosm_dataset$corridor_length,
                 colour = week_four_microcosm_dataset$temp_treatment)) +
  scale_shape_manual(values = c(3,4,15,17,19,16))+
  theme(legend.position = "none",
        aspect.ratio = 1)+
  labs(x = "PC1", y = "PC2")

#### week two PCA ####

#subset day 14 (halfway)
week_two_microcosm_dataset <- 
  microcosm_dataset %>% filter(day == 14)

#run the PCA
week_two_microcosm_PCA <- PCA(week_two_microcosm_dataset[,c(5:11)], graph = FALSE)

#examine eigenvalues to determine number of PCs
week_two_eig_vals <- get_eigenvalue(week_two_microcosm_PCA)
week_two_eig_vals #PCs 1, 2, 3, and 4 have eigenvalues greater than 1 and account for ~80% variance

#extract the results
week_two_vars <- get_pca_var(week_two_microcosm_PCA)

#correlation between variables and PCs
week_two_correlations <- week_two_vars$coord

#quality of representation
week_two_representation <- week_two_vars$cos2

#contribution of variables to PCs
week_two_contributions <- week_two_vars$contrib

#biplot of microcosms and species with confidence ellipses

week_two_biplot <- fviz_pca_biplot(week_two_microcosm_PCA,
                                    col.var = "black",
                                    label = "var",
                                    alpha.ind = 0.1,#make these points invisible to have diff shapes for points
                                    labelsize = 2,
                                    repel = TRUE,
                                    geom_ind = "points",
                                    col.ind = week_two_microcosm_dataset$temp_treatment,
                                    palette = c("#00AFBB", "#E7B800", "#FC4E07", "#0073C2FF"),
                                    mean.point = FALSE,
                                    addEllipses = TRUE,
                                    title = NULL,
                                    ellipse.type = "confidence",
                                    legend.title = "Temperature treatment",
                                    geom.ind = "point")+
  #add different shape points for long vs. short corridors (can't do this with fviz)
  geom_point(aes(shape = week_two_microcosm_dataset$corridor_length,
                 colour = week_two_microcosm_dataset$temp_treatment)) +
  scale_shape_manual(values = c(3,4,15,17,19,16))+
  theme(legend.position = "none",
        aspect.ratio = 1)+
  labs(x = "PC1", y = "PC2")


#### PERMANOVA and PERMDISP ####
library(vegan)

# week four 
#split the data into species and microcosm values
week_four_species_data <- week_four_microcosm_dataset[,c(5:10)]
week_four_microcosm_data <- week_four_microcosm_dataset[,c(1:2)]

#PERMANOVA 
adonis(week_four_species_data ~ corridor_length + temp_treatment,
       data = week_four_microcosm_data, method = "euclidian")

#post-hoc pairwise comparisons using pairwiseAdonis
library(devtools)
if (!require(pairwiseAdonis)) 
  install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

pairwise.adonis(week_four_species_data, week_four_microcosm_data$corridor_length) 
pairwise.adonis(week_four_species_data, week_four_microcosm_data$temp_treatment) 

#PERMDISP
#betadisper to test for differences in group homogeneities
#first calculate a distance matrix (done automatically in adonis)
week_four_distance_matrix <- vegdist(week_four_species_data, method = "bray")

#use the distance matrix to calculate multivariate dispersions
week_four_distance_model <- betadisper(week_four_distance_matrix, week_four_microcosm_data$corridor_length)

#use the dispersions to perform ANOVA showing if group dispersions are homogeneous
anova(week_four_distance_model)

#and again for temperature (can only do one factor at a time)
week_four_distance_model_temp <- betadisper(week_four_distance_matrix, week_four_microcosm_data$temp_treatment)
#use the dispersions to perform ANOVA showing if group dispersions are homogeneous
anova(week_four_distance_model_temp) #dispersions are heterogeneous, run Tukey to see which groups

TukeyHSD(week_four_distance_model_temp) #none of them (because of the adjustment)
#this plus the adonis shows that the groups differ in terms of their composition but not their dispersion
#meaning the compositions vary similarly 

# week two

#split the data into species and microcosm values

week_two_species_data <- week_two_microcosm_dataset[,c(5:11)]
week_two_microcosm_data <- week_two_microcosm_dataset[,c(1:2)]

#PERMANOVA
adonis(week_two_species_data ~ corridor_length + temp_treatment,
       data = week_two_microcosm_data, method = "euclidian")
pairwise.adonis(week_two_species_data, week_two_microcosm_data$corridor_length) 
pairwise.adonis(week_two_species_data, week_two_microcosm_data$temp_treatment)

#PERMDISP
#calculate the distance matrix
week_two_distance_matrix <- vegdist(week_two_species_data, method = "bray")

#calculate multivariate dispersions
week_two_distance_model_corridors <- betadisper(week_two_distance_matrix, week_two_microcosm_data$corridor_length)

#use the dispersions to perform ANOVA showing if group dispersions are homogeneous
anova(week_two_distance_model_corridors) #significant so the group dispersions are heterogeneous

#and again for temperature
week_two_distance_model_temp <- betadisper(week_two_distance_matrix, week_two_microcosm_data$temp_treatment)
#use the dispersions to perform ANOVA showing if group dispersions are homogeneous
anova(week_two_distance_model_temp) #dispersions are homogeneous
TukeyHSD(week_two_distance_model_temp)
#this plus the adonis shows that the groups differ in terms of their composition but not their dispersion
#meaning the compositions vary similarly 


