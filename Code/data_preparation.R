########################################################################################
### DATA LOADING AND PREPARATION ###
########################################################################################

### Preamble ###
library(openxlsx) # read.xlsx function
library(tidyverse) # data wrangling and plotting
library(vegan) # community analyses
library(betapart) # Beta diversity estimation
#---------------------------------------------------
#---------------------------------------------------

#set working directory and load data
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

data_raw <- openxlsx::read.xlsx("Data/ExpTempFluctuationsV2_Data_FINAL.xlsx", 
                                colNames = T, sheet = "Foglio1", detectDates = T) %>%
  mutate(across(Count_Blepharisma:Total_Bursaria,~replace(., is.na(.), 0))) %>% #replace count NAs with 0s
  mutate(Corridor = as.factor(Corridor), 
         Temp_Regime =  as.factor(Temp_Regime)) # ensure Corridor and Temp_Regime considered as factors

data_pooled <- data_raw %>% 
  dplyr::select(Date, Observer, Replicate, Temp_Regime, Corridor,Patch, Total_Blepharisma,
                Total_Didinium, Total_Colpidium, Total_Homolazoon, Total_Bursaria,
                Total_Spirostomum, Total_Paramecium) %>% # select variables of interest
  group_by(Date,Replicate,Temp_Regime,Corridor, Observer) %>% # group by categories for downstream processing
  dplyr::summarise_at(.vars = c("Total_Blepharisma", "Total_Didinium", "Total_Colpidium" ,
                                "Total_Homolazoon","Total_Bursaria", "Total_Spirostomum",
                                "Total_Paramecium"), .funs = c(sum ="sum")) %>% # sum associated patch abundances for each species
  tidyr::nest(data = Total_Blepharisma_sum:Total_Paramecium_sum) %>% # nest abundance columns to allow diversity metrics to be estimated for each row
  mutate(data = purrr::map(data, ~.x %>% # for each row:
                             mutate(Sh_div = vegan::diversity(.x), # estimate Shannon diversity
                                    Eveness = (Sh_div/vegan::specnumber(.x)), # species evenness
                                    Sp_rich = vegan::specnumber(.x)))) %>% # species richness
  tidyr::unnest(cols = c(data))%>% # unnest the dataframe 
  ungroup()%>%
  mutate(NumDays = as.numeric(difftime(Date, dplyr::first(Date), units = "days"))) # add time column representing "day of experiment"

data_patches<-data_raw%>%
  dplyr::select(Date, Replicate, Temp_Regime, Corridor,Patch, Total_Blepharisma,
                Total_Didinium, Total_Colpidium, Total_Homolazoon, Total_Bursaria,
                Total_Spirostomum, Total_Paramecium) %>% # select variables of interest
  group_by(Date,Replicate,Temp_Regime,Corridor) %>%
  tidyr::nest() %>% #nest abundance columns to allow diversity metrics to be estimated for each Patch combination
  mutate(data = purrr::map(data, ~.x %>% dplyr::select(-Patch))) %>% #drop patches
  mutate(data = purrr::map(data, ~.x %>%
                             mutate(beta = betapart::beta.multi.abund(.x)$beta.BRAY,# estimate Beta diversity
                                    Sp_rich = vegan::specnumber(.x),#species richness
                                    Patch = c("A","B"))))%>% # reintroduce patches
  tidyr::unnest(cols = c(data))%>% # unnest
  dplyr::ungroup() %>%
  mutate(NumDays = as.numeric(difftime(Date, dplyr::first(Date), units = "days"))) # add time column representing "day of experiment"

