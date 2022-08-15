#########################################################################################
### DISPERSAL TEST LONG VS SHORT CORRIDORS #####
#########################################################################################

#preable
library(openxlsx)
library(tidyverse) # data wrangling and plotting
library(ggplot2)
library(data.table)
library(lme4)


#load data
data<-read.xlsx("Data/Dispersal_experiments_data.xlsx", sheet = "Prey_species", colNames = T)

#perform Mann-Whitney-Wilcoxon Test
wilcox.test(N_individuals_2h ~ Corridor,data = data %>% filter(Species == "P.caudatum"))
wilcox.test(N_individuals_2h ~ Corridor,data = data %>% filter(Species == "B.japonicum"))
wilcox.test(N_individuals_2h ~ Corridor,data = data %>% filter(Species == "S.teres"))
wilcox.test(N_individuals_2h ~ Corridor,data = data %>% filter(Species == "Colpidium"))

#plot boxplots to see difference between long and short facet species
ggplot(data, aes (x = Corridor, y = N_individuals_2h, color = Species))+
  geom_boxplot()+
  facet_wrap(~Species, scales = "free")+
  theme_bw()+
  ggtitle("Prey species")
ggsave("dispersal_boxplots.pdf", units = "in", width = 10, height = 10)


###bootstrapping method
dataDT<-as.data.table(data)

data_pivot<-data %>% group_by(Species, Corridor) %>%
  summarise(avg = mean(N_individuals_2h)) %>% ungroup() %>% group_by(Species) %>% 
  pivot_wider(names_from = Corridor, values_from = avg) %>% 
  mutate(difference_in_average = Short - Long)
 

#create an empty dataset to fill with bootstrap permutations (using data.table)
difference.bootstrapped = data.table(Species = as.character(),
                                     N_individuals_short = as.numeric(),
                                     N_individuals_long = as.numeric(),
                                     diff_N_individuals = as.numeric())


#for loop to perform 1000 permutation bootstrapping: re sampling the data 1000 times and the calculate the mean of 
#dispersed individual in the simulated re sampling of the 5 replicates, for both long and short corridor
#and calculateadd a column where we calculate the difference between the means of short and long corridors 

for(i in 1:1000){
  dataDT_resampled = dataDT[
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
  difference.bootstrapped = rbind(difference.bootstrapped, difference.bootstrapped_current)
}
#calculate 95% confidence interval for the difference in the mean simulated dispersed individuals, for each species
quantile(difference.bootstrapped[Species == "P.caudatum"]$diff_N_individuals, c(0.025,0.975))
quantile(difference.bootstrapped[Species == "Colpidium"]$diff_N_individuals, c(0.025,0.975))
quantile(difference.bootstrapped[Species == "B.japonicum"]$diff_N_individuals, c(0.025,0.975))
quantile(difference.bootstrapped[Species == "S.teres"]$diff_N_individuals, c(0.025,0.975))



p_p.caudatum<-ggplot(difference.bootstrapped %>% filter(Species == "P.caudatum"),
                     aes(x = diff_N_individuals))+
  geom_histogram(colour="black", fill="white")+
  geom_vline(data = data_pivot %>% filter(Species == "P.caudatum"),
             aes(xintercept = difference_in_average),
             color="blue", linetype="dashed", size= 1)+
  geom_vline(aes(xintercept = -11.220),
             color = "red")+
  geom_vline(aes(xintercept = 49.205),
             color = "red")+
  labs(title = "P.caudatum")+
  xlim(-25,60)+
  theme_bw()

ggsave("p_p.caudatum.pdf", units = "in", width = 10, height = 10)

  
p_B.jap<-ggplot(difference.bootstrapped %>% filter(Species == "B.japonicum"),
                aes(x = diff_N_individuals))+
  geom_histogram(colour="black", fill="white")+
  geom_vline(data = data_pivot %>% filter(Species == "B.japonicum"),
             aes(xintercept = difference_in_average),
             color="blue", linetype="dashed", size= 1)+
  geom_vline(aes(xintercept =  16.2),
             color = "red")+
  geom_vline(aes(xintercept = 38.8),
             color = "red")+
  xlim(-25,60)+
  labs(title = "B.japonicum")+
  theme_bw()

ggsave("B.japonicum.pdf", units = "in", width = 10, height = 10)


p_Colp<-ggplot(difference.bootstrapped %>% filter(Species == "Colpidium"),
                aes(x = diff_N_individuals))+
  geom_histogram(colour="black", fill="white")+
  geom_vline(data = data_pivot %>% filter(Species == "Colpidium"),
             aes(xintercept = difference_in_average),
             color="blue", linetype="dashed", size= 1)+
  geom_vline(aes(xintercept =  4.0 ),
             color = "red")+
  geom_vline(aes(xintercept = 21.4),
             color = "red")+
  xlim(-25,60)+
  labs(title = "Colpidium")+
  theme_bw()

ggsave("Colpidium.pdf", units = "in", width = 10, height = 10)

  
p_Spir<-ggplot(difference.bootstrapped %>% filter(Species == "S.teres"),
               aes(x = diff_N_individuals))+
  geom_histogram(colour="black", fill="white")+
  geom_vline(data = data_pivot %>% filter(Species == "S.teres"),
             aes(xintercept = difference_in_average),
             color="blue", linetype="dashed", size= 1)+
  geom_vline(aes(xintercept =  7.795),
             color = "red")+
  geom_vline(aes(xintercept = 37.000),
             color = "red")+
  xlim(-25,60)+
  labs(title = "S.teres")+
  theme_bw()

ggsave("Spir.teres.pdf", units = "in", width = 10, height = 10)

  
  
  
#bootstrap
#create a datatable to contain the bootstrapped data and the difference between them
  
bootstrap.samples = data.table(bootstrap_id = as.integer(),
                               Species = as.character(),
                               Corridor = as.character(),
                               N_individuals_2h_resampled = as.integer())
for(i in 1:1000){
  bootstrap.samples_current = dataDT_resampled = dataDT[
    , list(bootstrap_id = i, N_individuals_2h_resampled = sample(x = N_individuals_2h, size = 5, replace = TRUE))
    , c("Species","Corridor")
  ]
  bootstrap.samples = rbind(bootstrap.samples, bootstrap.samples_current)
}




