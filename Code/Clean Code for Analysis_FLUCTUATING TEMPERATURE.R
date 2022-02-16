######Fluctuating temperature analysis#####
###Wolfe, Cerini, O'Brien, Besson, Clements 2022###

####load required packages####
library(tidyverse) 
library(viridis)
library(ggplot2)
library(ggpubr)
library(piecewiseSEM)
library(ggh4x)
library(rstatix)
library(openxlsx)
library(vegan)
library(entropart)
library(multcomp)
library(glmmTMB)
library(car)
library(fitdistrplus)
library(lmtest)
library(sjmisc)
library(splines)
#

########################################################################################
#### DATA LOADING AND PREPARATION ####
########################################################################################

###set working directory and load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data<-read.xlsx("ExpTempFluctuationsV2_Data_FINAL.xlsx", 
                colNames = T, sheet = "Foglio1", detectDates = T)
###put na as zeros
data[is.na(data)] <- 0

###group the data and select variables we are interested in
dataG<-data%>%
  group_by(Temp_Regime, Replicate, Date, Patch, Corridor) %>% 
  dplyr::select(Date, Observer, Replicate, Temp_Regime, Corridor,Patch, Total_Blepharisma,
                Total_Didinium, Total_Colpidium, Total_Homolazoon, Total_Bursaria,
                Total_Spirostomum, Total_Paramecium)


###combine data of the two patches 
data_sum_patch<-dataG %>% group_by(Date,Replicate,Temp_Regime,Corridor, Observer) %>%
  dplyr::summarise_at(.vars = c("Total_Blepharisma", "Total_Didinium", "Total_Colpidium" ,
                         "Total_Homolazoon","Total_Bursaria", "Total_Spirostomum",
                         "Total_Paramecium" ), .funs = c(sum ="sum"))

####calculating Shannon diversity index and add a measure of Eveness and Richness 
####and a column with number of days of the experiment
data_Sh_div<- data_sum_patch %>% 
  dplyr::ungroup() %>% 
  mutate(Sh_div = vegan::diversity(data_sum_patch[,6:12]))

data_Sh_div<-data_Sh_div %>% 
  mutate(Eveness = (Sh_div/specnumber(data_Sh_div[,6:12]))) %>% 
  mutate(Sp_rich = specnumber(data_sum_patch[,6:12])) %>% 
  mutate(NumDays = as.numeric(difftime(data_Sh_div$Date, data_Sh_div$Date[1], units = "days"))) 

###set these two columns as factors for needed in the following analysis
data_Sh_div$Corridor<-as.factor(data_Sh_div$Corridor)
data_Sh_div$Temp_Regime<-as.factor(data_Sh_div$Temp_Regime)



########################################################################################
##################### ANALYSIS AT THE END OF THE EXPERIMENT ############################
########################################################################################

#####################################################################################
#### ANALYSIS of SPECIES RICHNESS ####
#####################################################################################
##create a sub-dataset with just the last day of experiment
data_lastday<-data_Sh_div %>%
  filter(NumDays == 28) 

##set these two columns as factors for needed in the following analysis
data_lastday$Corridor<-as.factor(data_lastday$Corridor)
data_lastday$Temp_Regime<-as.factor(data_lastday$Temp_Regime)

##checking normality assumptions to run the two way anova
#Build the linear model
lm1<-lm(Sp_rich~ Corridor*Temp_Regime,
       data = data_lastday)

#Create a QQ plot of residuals
ggqqplot(residuals(lm1))

#Compute Shapiro-Wilk test of normality
shapiro_test(residuals(lm1))

#Check normality assumption by groups
data_lastday %>% 
  group_by(Temp_Regime, Corridor) %>%
  shapiro_test(Sp_rich)

#Create QQ plots for each cell of design
ggqqplot(data_lastday, "Sp_rich", ggtheme = theme_bw()) +
  facet_grid(Temp_Regime ~ Corridor)

#Homogneity of variance test
data_lastday %>% levene_test(Sp_rich~ Corridor*Temp_Regime)
####not all the subset of data are normally distributed but we are happy with the model and the other tests

##Run the anova (with two different function, don't know which one to keep, probably aov?)
summary(aov(Sp_rich~ Corridor*Temp_Regime,
            data = data_lastday))

##plotting the analysed data to see if the results match the data

ggplot(data_lastday, aes( x = Temp_Regime, y = Sp_rich, col = Temp_Regime))+
  geom_boxplot()+
  facet_grid(~Corridor)+
  ggtitle("Species richness last day")+
  theme_bw()





#################################################################################
#### ANALYSIS of Shannon Diversity ####
#################################################################################
##checking normality assumptions to run the two way anova
#Build the linear model#Create a QQ plot of residuals
lm2<-lm(Sh_div~ Corridor*Temp_Regime,
        data = data_lastday)
#Create a QQ plot of residuals
ggqqplot(residuals(lm2))

#Compute Shapiro-Wilk test of normality
shapiro_test(residuals(lm2))

#Check normality assumption by groups
data_lastday %>% 
  group_by(Temp_Regime, Corridor) %>%
  shapiro_test(Sh_div)
#Create QQ plots for each cell of design
ggqqplot(data_lastday, "Sh_div", ggtheme = theme_bw()) +
  facet_grid(Temp_Regime ~ Corridor)
##All the points fall approximately along the reference line, for each cell.
##Given all the other tests, we can assume normality of the data.

#Homogneity of variance testance
data_lastday %>% levene_test(Sh_div~ Corridor*Temp_Regime)



####Run the anova
anova_test(Sh_div~ Corridor*Temp_Regime,
           data = data_lastday)

####run multiple comparison post hoc test
data_lastday %>% pairwise_t_test(
  Sh_div ~ Temp_Regime, 
  p.adjust.method = "bonferroni"
)

##visualize the data 
#facet corridors
ggplot(data_lastday, aes( x = Temp_Regime, y = Sh_div, col = Temp_Regime))+
  geom_boxplot()+
  facet_grid(~Corridor)+
  ggtitle("Shannon diversity last day")+
  theme_bw()

##all together
ggplot(data_lastday, aes( x = Temp_Regime, y = Sh_div, col = Temp_Regime))+
  geom_boxplot()+
  ggtitle("Shannon diversity last day")+
  theme_bw()


#################################################################################
####ANALYSIS OF BETA DIVERSITY BETWEEN THE TWO PATCHES ####
###################################################################################
####Calculating Beta Diversity
##Create a new dataset with jsut the data we need 
data_Beta<-dataG%>%
  dplyr::select(Date, Replicate, Temp_Regime, Corridor,Patch, Total_Blepharisma,
                Total_Didinium, Total_Colpidium, Total_Homolazoon, Total_Bursaria,
                Total_Spirostomum, Total_Paramecium) %>% 
  dplyr::ungroup() %>% 
  group_by(Date,Replicate,Temp_Regime,Corridor)

##Calculate the Beta Diversity for each replica of each treatment

data_Beta <- data_Beta %>% group_by(Date,Replicate,Temp_Regime,Corridor) %>%
  tidyr::nest() %>%
  mutate(data = purrr::map(data, ~.x %>% dplyr::select(-Patch))) %>%
  mutate(data = purrr::map(data, ~.x %>%
                             mutate(beta = betapart::beta.multi.abund(.x)$beta.BRAY,
                                    Patch = c("A","B"))))%>%
  tidyr::unnest(cols = c(data))

##add a column with number of days 
data_Beta<-data_Beta %>% dplyr::ungroup() %>%
  mutate(NumDays = as.numeric(difftime(data_Beta$Date, data_Beta$Date[1], units = "days")))

##select just the last day to analyse the Beta diversity at the end of the experiment
dataBeta_lastday<-data_Beta %>% filter(Date =="2021-12-20") %>% filter((Patch=="A"))

#we need factors to avoid problems in the model
dataBeta_lastday$Temp_Regime<-as.factor(dataBeta_lastday$Temp_Regime)
dataBeta_lastday$Corridor<-as.factor(dataBeta_lastday$Corridor)

###checking normality assumptions to run the two way anova
#Build the linear model
dataBeta_lastday %>% ungroup()

modA<-lm(beta~ Corridor*Temp_Regime,
         data = dataBeta_lastday)

#Create a QQ plot of residuals
ggqqplot(residuals(modA))

#Compute Shapiro-Wilk test of normality
shapiro_test(residuals(modA))

#Check normality assumption by groups
dataBeta_lastday %>% 
  group_by(Temp_Regime, Corridor) %>%
  shapiro_test(beta)

#Create QQ plots for each cell of design
ggqqplot(dataBeta_lastday, "beta", ggtheme = theme_bw()) +
  facet_grid(Temp_Regime ~ Corridor)

#Homogneity of variance test
dataBeta_lastday %>% ungroup() %>%  
  levene_test(beta ~ Corridor*Temp_Regime)

####Run the anova
summary(aov(modA))

####run multiple comparison post hoc test
pairwise_t_test(dataBeta_lastday, beta ~ Temp_Regime, 
                p.adjust.method = "bonferroni")

##visualize data
ggplot(dataBeta_lastday, aes( x = Temp_Regime, y = beta, col = Temp_Regime))+
  geom_boxplot()+
  facet_grid(~Corridor)+
  ggtitle("Beta diversity between two patches LAST DAY")+
  theme_bw()

ggplot(dataBeta_lastday, aes( x = Temp_Regime, y = beta, col = Temp_Regime))+
  geom_boxplot()+
  ggtitle("Beta diversity between two patches LAST DAY")+
  theme_bw()








########################################################################################
##########################ANALYSIS THROUGH TIME#########################################
########################################################################################

################################################################################
#### ANALYSIS of SPECIES RICHNESS THROUHG TIME#### 
################################################################################

####Model Building

##model with all interactions with the time factor (in this case NumDays:number of days of the experiment)
mod_Sp_time<-glm(Sp_rich~NumDays +
                   Temp_Regime + Corridor +NumDays:Temp_Regime + NumDays:Corridor + Corridor:Temp_Regime+
                   Temp_Regime:Corridor:NumDays,
                 data = data_Sh_div,
                 family = "poisson" )


##visual checking of the model
plot(mod_Sp_time)

##plot the residual of the model to look at their distribution
hist(resid(mod_Sp_time))

##check normality of the residuals with Shapiro-Wilk test
shapiro_test(resid(mod_Sp_time))

##model outputs
summary(mod_Sp_time)
Anova(mod_Sp_time)


##model simplification 1 (without cubic interaction)
mod_simp1<-glm(Sp_rich~NumDays +
                 Temp_Regime + Corridor +NumDays:Temp_Regime + NumDays:Corridor + Corridor:Temp_Regime,
               data = data_Sh_div,
              family = "poisson" )


##visual checking of the model
plot(mod_simp1)

###plot the residual of the model to look at their distribution
hist(resid(mod_simp1))

##check normality of the residuals with Shapiro-Wilk test
shapiro.test(resid(mod_simp1))


##model outputs
summary(mod_simp1)
Anova(mod_simp1)


##model simplification 2 (without quadratic interactions)
mod_simp2<-glm(Sp_rich~NumDays +
                 Temp_Regime + Corridor,
               data = data_Sh_div,
               family = "poisson" )


##visual checking of the model
plot(mod_simp2)

###plot the residual of the model to look at their distribution

hist(resid(mod_simp2))

##check normality of the residuals with Shapiro-Wilk test
shapiro.test(resid(mod_simp2))


##model output
summary(mod_simp2)
Anova(mod_simp2)

##model simplification 3 (without "corridor" factor)
mod_simp3<-glm(Sp_rich~NumDays +
                 Temp_Regime,
               data = data_Sh_div,
               family = "poisson" )

##control model just with time to be sure of the significant effect of the Temperature Regime Factor
mod_simp4<-glm(Sp_rich~NumDays,
               data = data_Sh_div,
               family = "poisson" )

##visual checking of the model
plot(mod_simp3)


###plot the residual of the model to look at their distribution
hist(resid(mod_simp3))

##check normality of the residuals with Shapiro-Wilk test
shapiro.test(resid(mod_simp3))


##model output
summary(mod_simp3)
Anova(mod_simp3)

##multiple comparison post hoc analysis
summary(glht(mod_simp3, mcp(Temp_Regime = "Tukey")))

##Likelihood Ratio Test to compare if our model selection is significant 
lrtest(mod_Sp_time,mod_simp1)
lrtest(mod_simp1,mod_simp2)
lrtest(mod_simp2,mod_simp3)
lrtest(mod_simp3,mod_simp4)




##predict values from the model and plot it, to see the fit of the model to real data
df.simp<-data.frame(preds = predict(mod_simp3, newdata = data_Sh_div,type = "response" ),
                    upr = predict(mod_simp3, newdata = data_Sh_div ,se.fit = T,type = "response" )$fit + predict(mod_simp3, newdata = data_Sh_div,se.fit = T,type = "response" )$se.fit ,
                    lwr = predict(mod_simp3, newdata = data_Sh_div ,se.fit = T ,type = "response")$fit - predict(mod_simp3, newdata = data_Sh_div ,se.fit = T ,type = "response")$se.fit , data_Sh_div)

ggplot(df.simp, aes(x =NumDays, y= Sp_rich,
                            col =Temp_Regime, fill = Temp_Regime))+
  geom_point(position=position_jitter(width=0.3, height=0.05, seed = 1), alpha = .2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),col="transparent",alpha=0.1)+
  #facet_wrap(~Corridor)+
  ggtitle("Species Richness Poisson model")+
  theme_bw()




################################################################################
#### ANALYSIS of SHANNON DIVERSITY THROUHG TIME ####
################################################################################
##plot the data to see distribution
data_Sh_div %>% group_by(Replicate, Corridor, Temp_Regime)

####Build the the model with cubic Spline for the time variable 
##model with all interactions with the time factor (in this case NumDays:number of days of the experiment)
mod_SH_time1_SPLI<-glm(Sh_div~ ns(NumDays,3) +
                         Temp_Regime + 
                         Corridor + 
                         ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor + Temp_Regime:Corridor+
                         ns(NumDays,3):Temp_Regime:Corridor, data =data_Sh_div, family = "gaussian")

##visual checking of the model
plot(mod_SH_time1_SPLI)   

##plot the residual of the model to look at their distribution
hist(residuals(mod_SH_time1_SPLI))

##check normality of the residuals with Shapiro-Wilk test
shapiro.test(residuals(mod_SH_time1_SPLI))

##Output of the model
summary(mod_SH_time1_SPLI)
Anova(mod_SH_time1_SPLI, test.statistic = "F")


####no interaction of the cubic term, we take it out of the model and re run
##model simplification 1 (deleting cubic term)
mod_SH_time1_SPLI_Simp1<-glm(Sh_div~ ns(NumDays,3) +
                               Temp_Regime + 
                               Corridor + 
                               ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor + Temp_Regime:Corridor,
                             data =data_Sh_div, family = "gaussian")

##visual checking of the model
plot(mod_SH_time1_SPLI_Simp1)   

##plot the residual of the model to look at their distribution
hist(residuals(mod_SH_time1_SPLI_Simp1))

##check normality of the residuals with Shapiro-Wilk test
shapiro.test(residuals(mod_SH_time1_SPLI_Simp1))

###output of the model
summary(mod_SH_time1_SPLI_Simp1)
Anova(mod_SH_time1_SPLI_Simp1, test.statistic = "F")


####no effect of the Temp_Regime:Corridor interaction, we take it out of the model and re run
##model simplification 2 
mod_SH_time1_SPLI_Simp2<-glm(Sh_div~ ns(NumDays,3) +
                               Temp_Regime + 
                               Corridor + 
                               ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor,
                             data =data_Sh_div, family = "gaussian")


##visual checking of the model
plot(mod_SH_time1_SPLI_Simp2)   

##plot the residual of the model to look at their distribution
hist(residuals(mod_SH_time1_SPLI_Simp2))

##check normality of the residuals with Shapiro-Wilk test
shapiro.test(residuals(mod_SH_time1_SPLI_Simp2))

###results of the model
summary(mod_SH_time1_SPLI_Simp2)
Anova(mod_SH_time1_SPLI_Simp2, test.statistic = "F")

##plot the average difference in Shannon diversity across treatments to see why we have a significant
##role of temp_regime

ggplot(data_Sh_div, aes(x = Temp_Regime, y = Sh_div, col = Temp_Regime))+
  geom_boxplot()+
  ggtitle("Shannon diversity all time points together")+
  theme_bw()

###visually it does not look like we have on average big difference as the model ANOVA is telling
##running an anova to check confirm that there is no significant difference
aov_mod_SH_time1_SPLI_Simp2<-aov(Sh_div~
                 Temp_Regime,
               data =data_Sh_div)
summary(aov_mod_SH_time1_SPLI_Simp2)



##Likelihood Ratio Test to compare if our model selection is significant 
lrtest(mod_SH_time1_SPLI_noint,mod_SH_time1_SPLI_Simp1)
lrtest(mod_SH_time1_SPLI_Simp1,mod_SH_time1_SPLI_Simp2)



##predict values from the model and plot it, to see the fit of the model to real data
gg.df_noint <- data.frame(preds = predict(mod_SH_time1_SPLI_Simp2, newdata = data_Sh_div,type = "response" ),
                    upr = predict(mod_SH_time1_SPLI_Simp2, newdata = data_Sh_div ,se.fit = T,type = "response" )$fit +
                      predict(mod_SH_time1_SPLI_Simp2, newdata = data_Sh_div,se.fit = T,type = "response" )$se.fit ,
                    lwr = predict(mod_SH_time1_SPLI_Simp2, newdata = data_Sh_div ,se.fit = T ,type = "response")$fit -
                      predict(mod_SH_time1_SPLI_Simp2, newdata = data_Sh_div ,se.fit = T ,type = "response")$se.fit , data_Sh_div)


ggplot(gg.df_noint, aes(x =NumDays, y= Sh_div,group = Temp_Regime, col =Temp_Regime, fill = Temp_Regime))+
  geom_point(alpha = 0.2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  facet_grid(~Corridor)+
  theme_bw()

ggplot(gg.df_noint, aes(x =NumDays, y= Sh_div,col =Temp_Regime, fill = Temp_Regime,
                        linetype = Corridor))+
  geom_point(alpha = 0.2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  theme_bw()

###fit a model just with the significant interaction of TEMP TEGIME with time to 
##plot the differences
modTemp<-glm(Sh_div~ ns(NumDays,3)+Temp_Regime + 
              + ns(NumDays,3):Temp_Regime,
            data =data_Sh_div, family = "gaussian")

gg.df_TEMP <- data.frame(preds = predict(modTemp, newdata = data_Sh_div,type = "response" ),
                         upr = predict(modTemp, newdata = data_Sh_div ,se.fit = T,type = "response" )$fit +
                           predict(modTemp, newdata = data_Sh_div,se.fit = T,type = "response" )$se.fit ,
                         lwr = predict(modTemp, newdata = data_Sh_div ,se.fit = T ,type = "response")$fit -
                           predict(modTemp, newdata = data_Sh_div ,se.fit = T ,type = "response")$se.fit , data_Sh_div)


ggplot(gg.df_TEMP, aes(x =NumDays, y= Sh_div,group = Temp_Regime,
                       col =Temp_Regime, fill = Temp_Regime))+
  geom_point(alpha = 0.2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  theme_bw()



################################################################################
#### ANALYSIS of Shannon diversity VARIANCE THROUHG TIME ####
################################################################################
###calculating median Sh_diversity between replicates + upper and lower standard deviations
###+ variance

data_Sh_div_VAR<-data_Sh_div %>%
  dplyr::group_by(Date,Temp_Regime, Corridor, NumDays) %>% 
  dplyr::summarise(median.Sh_div = median(Sh_div), sd.val = sd(Sh_div), Sh_variance = var(Sh_div),
                   upr = median.Sh_div + sd.val, lwr = median.Sh_div - sd.val)


####avoid considering the first day in which we have variance zero, we create a dataframe without initial day
data_Sh_div_VAR_no_ZERO<-data_Sh_div_VAR %>%filter(NumDays != 0)

##put the two variables as factors
data_Sh_div_VAR_no_ZERO$Temp_Regime<-as.factor(data_Sh_div_VAR_no_ZERO$Temp_Regime)
data_Sh_div_VAR_no_ZERO$Corridor<-as.factor(data_Sh_div_VAR_no_ZERO$Corridor)

####Build the the model with cubic Spline for the time variable 
##model with all interactions with the time factor (in this case NumDays:number of days of the experiment)
mod_SH_VAR_SPL1<-glm(Sh_variance~ ns(NumDays,3) +
                               Temp_Regime + 
                               Corridor + 
                               ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor + Temp_Regime:Corridor+
                               ns(NumDays,3):Temp_Regime:Corridor, data =data_Sh_div_VAR_no_ZERO,
                     family = gaussian(link = "identity"))



##visual checking of the model
plot(mod_SH_VAR_SPL1)   

##plot the residual of the model to look at their distribution
hist(residuals(mod_SH_VAR_SPL1))

##check normality of the residuals with Shapiro-Wilk test
shapiro.test(residuals(mod_SH_VAR_SPL1))

###Output of the model

summary(mod_SH_VAR_SPL1)
Anova(mod_SH_VAR_SPL1, test.statistic = "F")


##post hoc analysis for interaction factor Temp_Regime:Corridor
anova_var<-aov(Sh_variance~
      Temp_Regime*Corridor,
      data =data_Sh_div_VAR_no_ZERO)
summary(anova_var)

###extract pairwise comparison and create a dataframe with just the significant p values of comparisons
sig_comp<-as.data.frame(TukeyHSD(anova_var)$"Temp_Regime:Corridor")
sig_comp<-sig_comp %>% filter(`p adj` < 0.05)
sig_comp

##predict values from the model and plot it, to see the fit of the model to real data
gg.df_VAR <- data.frame(preds = predict(mod_SH_VAR_SPL1, newdata = data_Sh_div_VAR_no_ZERO,type = "response" ),
                          upr = predict(mod_SH_VAR_SPL1, newdata = data_Sh_div_VAR_no_ZERO ,se.fit = T,type = "response" )$fit +
                            predict(mod_SH_VAR_SPL1, newdata = data_Sh_div_VAR_no_ZERO,se.fit = T,type = "response" )$se.fit ,
                          lwr = predict(mod_SH_VAR_SPL1, newdata = data_Sh_div_VAR_no_ZERO ,se.fit = T ,type = "response")$fit -
                            predict(mod_SH_VAR_SPL1, newdata = data_Sh_div_VAR_no_ZERO ,se.fit = T ,type = "response")$se.fit , data_Sh_div_VAR_no_ZERO)


ggplot(gg.df_VAR, aes(x =NumDays, y= Sh_variance,group = Temp_Regime, col =Temp_Regime, fill = Temp_Regime))+
  geom_point(alpha = 0.2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  facet_grid(~Corridor)+
  theme_bw()

ggplot(data_Sh_div_VAR_no_ZERO, aes(x = Temp_Regime, y = Sh_variance, col = Corridor))+
  geom_boxplot()+
  theme_bw()


##code the simpler model with splines
mod_SH_VAR_SPL2<-glm(Sh_variance~ ns(NumDays,3) +
                       Temp_Regime + 
                       Corridor + 
                       ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor + Temp_Regime:Corridor,
                     data =data_Sh_div_VAR_no_ZERO,
                     family = gaussian(link = "identity"))


##Likelihood Ratio Test to compare if our model selection is significant 
lrtest(mod_SH_VAR_SPL2,mod_SH_VAR_SPL1)


#####end####