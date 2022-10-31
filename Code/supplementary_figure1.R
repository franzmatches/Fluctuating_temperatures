#########################################################################################
### DISPERSAL TEST LONG VS SHORT CORRIDORS #####
#########################################################################################

#preamble
library(openxlsx)
library(tidyverse) 
library(ggplot2)
library(data.table)
library(ggpubr)

####PREY SPECIES####
#load data
dispersal_data_prey<-read.xlsx("Data/Dispersal_experiments_data.xlsx", sheet = "Prey_species", colNames = T)

#pivot data to display the mean value for each species and corridor length,
#and the value of the difference between the two
dispersal_data_prey_pivot <-dispersal_data_prey %>% group_by(Species, Corridor) %>%
  summarise(avg = mean(N_individuals_2h)) %>% ungroup() %>% group_by(Species) %>% 
  pivot_wider(names_from = Corridor, values_from = avg) %>% 
  mutate(difference_in_average = Short - Long)


#create an empty dataset to fill with bootstrap permutations (using data.table)
prey.difference.bootstrapped <- data.table::data.table(Species = as.character(),
                                     N_individuals_short = as.numeric(),
                                     N_individuals_long = as.numeric(),
                                     diff_N_individuals = as.numeric())

#for loop to perform 1000 permutation bootstrapping: re sampling the data 1000 times and then calculate the mean of 
#dispersed individual in the simulated re sampling of the 5 replicates, for both long and short corridor,
#and add a column where we calculate the difference between the means of short and long corridors 

#put the data in Data.table format
dispersal_prey_dataDT<-as.data.table(dispersal_data_prey)
for(i in 1:1000){
  dataDT_resampled = dispersal_prey_dataDT[
    , list(N_individuals_2h_resampled = sample(x = N_individuals_2h, size = 5, replace = TRUE))
    , c("Species","Corridor")
  ]
  
  difference.bootstrapped_current = merge(dataDT_resampled[Corridor == "Short",
                                                           list(N_individuals_short = mean(N_individuals_2h_resampled)),
                                                           c("Species")],
                                          dataDT_resampled[Corridor == "Long",
                                                           list(N_individuals_long = mean(N_individuals_2h_resampled)),
                                                           c("Species")],
                                          by = "Species")
  
  difference.bootstrapped_current[, diff_N_individuals := N_individuals_short - N_individuals_long]
  prey.difference.bootstrapped = rbind(prey.difference.bootstrapped, difference.bootstrapped_current)
}

#calculate 95% confidence interval for the difference in the mean simulated dispersed individuals, for each species
p.caudatum_quantile <- quantile(prey.difference.bootstrapped[Species == "P.caudatum"]$diff_N_individuals, c(0.025,0.975))
colpidium_quantile <- quantile(prey.difference.bootstrapped[Species == "Colpidium"]$diff_N_individuals, c(0.025,0.975))
b.japonicum_quantile <- quantile(prey.difference.bootstrapped[Species == "B.japonicum"]$diff_N_individuals, c(0.025,0.975))
s.teres_quantile <- quantile(prey.difference.bootstrapped[Species == "S.teres"]$diff_N_individuals, c(0.025,0.975))


p_p.caudatum<-ggplot(prey.difference.bootstrapped %>% filter(Species == "P.caudatum"),
                     aes(x = diff_N_individuals))+
  geom_histogram(colour="black", fill="white")+
  geom_vline(data = dispersal_data_prey_pivot %>% filter(Species == "P.caudatum"),
             aes(xintercept = difference_in_average),
             color="blue", linetype="dashed", size= 1)+
  geom_vline(aes(xintercept = p.caudatum_quantile[1]),
             color = "red")+
  geom_vline(aes(xintercept = p.caudatum_quantile[2]),
             color = "red")+
  labs(title = "P.caudatum")+
  xlim(-25,60)+
  theme_bw(base_size = 20)

p_B.jap<-ggplot(prey.difference.bootstrapped %>% filter(Species == "B.japonicum"),
                aes(x = diff_N_individuals))+
  geom_histogram(colour="black", fill="white")+
  geom_vline(data = dispersal_data_prey_pivot %>% filter(Species == "B.japonicum"),
             aes(xintercept = difference_in_average),
             color="blue", linetype="dashed", size= 1)+
  geom_vline(aes(xintercept =  b.japonicum_quantile[1]),
             color = "red")+
  geom_vline(aes(xintercept = b.japonicum_quantile[2]),
             color = "red")+
  xlim(-25,60)+
  labs(title = "B.japonicum")+
  theme_bw(base_size = 20)


p_Colp<-ggplot(prey.difference.bootstrapped %>% filter(Species == "Colpidium"),
               aes(x = diff_N_individuals))+
  geom_histogram(colour="black", fill="white")+
  geom_vline(data = dispersal_data_prey_pivot %>% filter(Species == "Colpidium"),
             aes(xintercept = difference_in_average),
             color="blue", linetype="dashed", size= 1)+
  geom_vline(aes(xintercept =  colpidium_quantile[1] ),
             color = "red")+
  geom_vline(aes(xintercept = colpidium_quantile[2]),
             color = "red")+
  xlim(-25,60)+
  labs(title = "C.striatum")+
  theme_bw(base_size = 20)


p_Spir<-ggplot(prey.difference.bootstrapped %>% filter(Species == "S.teres"),
               aes(x = diff_N_individuals))+
  geom_histogram(colour="black", fill="white")+
  geom_vline(data = dispersal_data_prey_pivot %>% filter(Species == "S.teres"),
             aes(xintercept = difference_in_average),
             color="blue", linetype="dashed", size= 1)+
  geom_vline(aes(xintercept =  s.teres_quantile[1]),
             color = "red")+
  geom_vline(aes(xintercept = s.teres_quantile[2]),
             color = "red")+
  xlim(-25,60)+
  labs(title = "S.teres")+
  theme_bw(base_size = 20)

####PREDATOR SPECIES####
#load data
dispersal_data_predator<-read.xlsx("Data/Dispersal_experiments_data.xlsx", sheet = "Predator_species", colNames = T)

#pivot data to display the mean value for each species and corridor length,
#and the value of the difference between the two
dispersal_data_predator_pivot <-dispersal_data_predator %>% group_by(Species, Corridor) %>%
  summarise(avg = mean(N_individuals_4h)) %>% ungroup() %>% group_by(Species) %>% 
  pivot_wider(names_from = Corridor, values_from = avg) %>% 
  mutate(difference_in_average = Short - Long)

#create an empty dataset to fill with bootstrap permutations (using data.table)
predator.difference.bootstrapped = data.table(Species = as.character(),
                                     N_individuals_short = as.numeric(),
                                     N_individuals_long = as.numeric(),
                                     diff_N_individuals = as.numeric())
#for loop to perform 1000 permutation bootstrapping: re sampling the data 1000 times and then calculate the mean of 
#dispersed individual in the simulated re sampling of the 5 replicates, for both long and short corridor,
#and add a column where we calculate the difference between the means of short and long corridors  

#put the data in Data.table format
dispersal_predator_dataDT<-as.data.table(dispersal_data_predator)
for(i in 1:1000){
  dataDT_resampled = dispersal_predator_dataDT[
    , list(N_individuals_4h_resampled = sample(x = N_individuals_4h, size = 5, replace = TRUE))
    , c("Species","Corridor")
  ]
  
  difference.bootstrapped_current = merge(dataDT_resampled[Corridor == "Short",
                                                           list(N_individuals_short = mean(N_individuals_4h_resampled)),
                                                           c("Species")],
                                          dataDT_resampled[Corridor == "Long",
                                                           list(N_individuals_long = mean(N_individuals_4h_resampled)),
                                                           c("Species")],
                                          by = "Species")
  
  difference.bootstrapped_current[, diff_N_individuals := N_individuals_short - N_individuals_long]
  predator.difference.bootstrapped = rbind(predator.difference.bootstrapped, difference.bootstrapped_current)
}
#calculate 95% confidence interval for the difference in the mean simulated dispersed individuals
d.nasutum_quantile <- quantile(predator.difference.bootstrapped[Species == "Didinium_nasutum"]$diff_N_individuals, c(0.025,0.975))

#plot histograms of the bootstrapping simulated differences between long and short corridor 
p_Didinium_nasutum<-ggplot(predator.difference.bootstrapped %>% filter(Species == "Didinium_nasutum"),
                           aes(x = diff_N_individuals))+
  geom_histogram(colour="black", fill="white")+
  geom_vline(data = dispersal_data_predator_pivot %>% filter(Species == "Didinium_nasutum"),
             aes(xintercept = difference_in_average),
             color="blue", linetype="dashed", size= 1)+
  geom_vline(aes(xintercept = d.nasutum_quantile[1]),
             color = "red")+
  geom_vline(aes(xintercept = d.nasutum_quantile[2]),
             color = "red")+
  labs(title = "D.nasutum")+
  theme_bw(base_size = 20)

####Merge Plots####

dispersal_plot <- ggpubr::ggarrange(
  ggpubr::ggarrange(p_B.jap,p_Colp,p_Spir,p_p.caudatum,ncol = 2,nrow=2),
  ggpubr::ggarrange(NULL, p_Didinium_nasutum, NULL, ncol = 3, widths = c(1,2,1)),
  ncol = 1,heights = c(2,1))

ggsave("Results/Supplementary_Figure1.tiff", dispersal_plot, units = "in", width = 10, height = 10)



