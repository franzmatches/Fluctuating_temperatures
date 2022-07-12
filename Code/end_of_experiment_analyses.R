#########################################################################################
### ANALYSIS AT THE END OF THE EXPERIMENT ###
#########################################################################################

### Preamble ###
library(tidyverse) # data wrangling and plotting
library(ggpubr) # multi panel figures
library(rstatix) # pipe-friendly model assumption tests
library(multcomp) # parametric model comparison functions
library(lmtest) # model comparisons

source("Code/data_preparation.R")
palette <- c("#440154FF", "#3B528BFF","#21908CFF", "#5DC863FF")

#create a sub-dataset with just the last day of experiment
data_pooled_lastday <- data_pooled %>%
  filter(NumDays ==max(NumDays))

#select just the last day to analyse the Beta diversity at the end of the experiment
data_patches_lastday <- data_patches %>% 
  filter(NumDays ==max(NumDays) & Patch == "A")

#----------------------------------------------------------------------------------------
#### ANALYSIS of Shannon Diversity ####
#----------------------------------------------------------------------------------------

##checking normality assumptions to run the two way anova
#Build the linear model#Create a QQ plot of residuals
lm2<-lm(Sh_div~ Corridor*Temp_Regime,
        data = data_pooled_lastday)
#Create a QQ plot of residuals
ggpubr::ggqqplot(residuals(lm2))

#Compute Shapiro-Wilk test of normality
rstatix::shapiro_test(residuals(lm2))

#Check normality assumption by groups
data_pooled_lastday %>% 
  group_by(Temp_Regime, Corridor) %>%
  rstatix::shapiro_test(Sh_div)

#Create QQ plots for each cell of design
ggpubr::ggqqplot(data_pooled_lastday, "Sh_div", ggtheme = theme_bw()) +
  facet_grid(Temp_Regime ~ Corridor)
##All the points fall approximately along the reference line, for each call.
##Given all the other tests, we can assume normality of the data.

#Homogneity of variance testance
data_pooled_lastday %>% rstatix::levene_test(Sh_div~ Corridor*Temp_Regime)

####Run the anova
rstatix::anova_test(Sh_div~ Corridor*Temp_Regime,
                    data = data_pooled_lastday)

####run multiple comparison post hoc test
data_pooled_lastday %>% rstatix::pairwise_t_test(
  Sh_div ~ Temp_Regime, 
  p.adjust.method = "bonferroni")

##visualize the data by corridor
ggplot(data_pooled_lastday, aes( x = Temp_Regime, y = Sh_div, col = Temp_Regime))+
  geom_boxplot()+
  facet_grid(~Corridor)+
  ggtitle("Shannon diversity last day")+
  theme_bw()

##all together
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

#----------------------------------------------------------------------------------------
## ANALYSIS OF BETA DIVERSITY BETWEEN THE TWO PATCHES ##
#----------------------------------------------------------------------------------------

###checking normality assumptions to run the two way anova
#Build the linear model
modA<-lm(beta~ Corridor*Temp_Regime,
         data = data_patches_lastday)

#Create a QQ plot of residuals
ggpubr::ggqqplot(residuals(modA))

#Compute Shapiro-Wilk test of normality
rstatix::shapiro_test(residuals(modA))

#Check normality assumption by groups
data_patches_lastday %>% 
  group_by(Temp_Regime, Corridor) %>%
  rstatix::shapiro_test(beta)

#Create QQ plots for each cell of design
ggpubr::ggqqplot(data_patches_lastday, "beta", ggtheme = theme_bw()) +
  facet_grid(Temp_Regime ~ Corridor)

#Homogneity of variance test
data_patches_lastday %>% ungroup() %>%  
  rstatix::levene_test(beta ~ Corridor*Temp_Regime)

####Run the anova
summary(aov(modA))

####run multiple comparison post hoc test
rstatix::pairwise_t_test(data_patches_lastday, beta ~ Temp_Regime, 
                         p.adjust.method = "bonferroni")

##visualize data
ggplot(data_patches_lastday, aes( x = Temp_Regime, y = beta, col = Temp_Regime))+
  geom_boxplot()+
  facet_grid(~Corridor)+
  ggtitle("Beta diversity between two patches LAST DAY")+
  theme_bw()

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

end_plots <- ggarrange(shannon_end, beta_end, align = "hv", labels = "auto", label.x = 0.08, label.y = 0.71)

ggsave("end_plots.tiff", end_plots)

#----------------------------------------------------------------------------------------
## ANALYSIS of SPECIES RICHNESS ##
#----------------------------------------------------------------------------------------

##checking normality assumptions
#Build the linear model
lm1<-lm(Sp_rich~ Corridor*Temp_Regime,
        data = data_pooled_lastday)

#Create a QQ plot of residuals
ggpubr::ggqqplot(residuals(lm1))

#Compute Shapiro-Wilk test of normality
shapiro.test(residuals(lm1))

#Check normality assumption by groups
data_pooled_lastday %>% 
  group_by(Temp_Regime, Corridor) %>%
  rstatix::shapiro_test(Sp_rich)

#Create QQ plots for each cell of design
ggpubr::ggqqplot(data_pooled_lastday, "Sp_rich", ggtheme = theme_bw()) +
  facet_grid(Temp_Regime ~ Corridor)

#Homogneity of variance test
data_pooled_lastday %>% rstatix::levene_test(Sp_rich~ Corridor*Temp_Regime)
#not all the subset of data are normally distributed but we are happy with the model and the other tests

##Run the anova (with two different function, don't know which one to keep, probably aov?)
summary(aov(Sp_rich~ Corridor*Temp_Regime,
            data = data_pooled_lastday))

##plotting the analysed data to see if estimates match the data

ggplot(data_pooled_lastday, aes( x = Temp_Regime, y = Sp_rich, col = Temp_Regime))+
  geom_boxplot()+
  facet_grid(~Corridor)+
  ggtitle("Species richness last day")+
  theme_bw()

