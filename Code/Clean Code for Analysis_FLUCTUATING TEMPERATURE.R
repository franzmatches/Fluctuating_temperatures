######Fluctuating temperature analysis####
###Wolfe, Cerini, O'Brien, Besson, Clements 2022###
rm(list=ls(all=TRUE))

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
library(DHARMa)
#####


###set working directory and load data
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data<-read.xlsx("ExpTempFluctuationsV2_Data_FINAL.xlsx", 
                colNames = T, sheet = "Foglio1", detectDates = T)
###put na as zeros
data[is.na(data)] <- 0
dataG<-data%>%
  group_by(Temp_Regime, Replicate, Date, Patch, Corridor) %>% 
  dplyr::select(Date, Observer, Replicate, Temp_Regime, Corridor,Patch, Total_Blepharisma,
                Total_Didinium, Total_Colpidium, Total_Homolazoon, Total_Bursaria,
                Total_Spirostomum, Total_Paramecium)


###combine data of the two patches 
colnames((dataG))
data_sum_patch<-dataG %>% group_by(Date,Replicate,Temp_Regime,Corridor, Observer) %>%
  dplyr::summarise_at(.vars = c("Total_Blepharisma", "Total_Didinium", "Total_Colpidium" ,
                         "Total_Homolazoon","Total_Bursaria", "Total_Spirostomum",
                         "Total_Paramecium" ), .funs = c(sum ="sum"))

####calculating Shannon diversity index and add a measure of Eveness and Richness and a colummn with number of days
data_Sh_div<- data_sum_patch %>% 
  dplyr::ungroup() %>% 
  mutate(Sh_div = vegan::diversity(data_sum_patch[,6:12]))

data_Sh_div<-data_Sh_div %>% 
  mutate(Eveness = (Sh_div/specnumber(data_Sh_div[,6:12]))) %>% 
  mutate(Sp_rich = specnumber(data_sum_patch[,6:12])) %>% 
  mutate(NumDays = as.numeric(difftime(data_Sh_div$Date, data_Sh_div$Date[1], units = "days"))) 


data_Sh_div$Corridor<-as.factor(data_Sh_div$Corridor)
data_Sh_div$Temp_Regime<-as.factor(data_Sh_div$Temp_Regime)


#ANALYSIS of SPECIES RICHNESS LAST DAY WITH INTERACTION Corridor*Temp_Regime##########
###create a subdataset with just the last day of experiment
data_lastday<-data_Sh_div %>%
  filter(NumDays == 28) 

data_lastday$Corridor<-as.factor(data_lastday$Corridor)
data_lastday$Temp_Regime<-as.factor(data_lastday$Temp_Regime)

###checking assumptions to run the two way anova
lm<-lm(Sp_rich~ Corridor*Temp_Regime,
       data = data_lastday)
ggqqplot(residuals(lm))

##normality of data and model
data_lastday %>% 
  group_by(Temp_Regime, Corridor) %>%
  shapiro_test(Sp_rich)


ggqqplot(data_lastday, "Sp_rich", ggtheme = theme_bw()) +
  facet_grid(Temp_Regime ~ Corridor)

###Homogneity of variance test
data_lastday %>% levene_test(Sp_rich~ Corridor*Temp_Regime)
####not all the subset of data are normally distributed but we are happy with the model

##run the anova (with two different function, don't know chich one to keep, probably aov?)
summary(aov(Sp_rich~ Corridor*Temp_Regime,
            data = data_lastday))

anova_test(Sp_rich~ Corridor*Temp_Regime,
    data = data_lastday)


####plotting the analysed data to see if the results makes sense

ggplot(data_lastday, aes( x = Temp_Regime, y = Sp_rich, col = Temp_Regime))+
  geom_boxplot()+
  facet_grid(~Corridor)+
  ggtitle("Species richness last day")+
  theme_bw()



##### Shannon Diversity Index THE LAST DAY ACROSS THE TREATMENTS######

##visualize the data to see what to expect
##with corridors
ggplot(data_lastday, aes( x = Temp_Regime, y = Sh_div, col = Temp_Regime))+
  geom_boxplot()+
  facet_grid(~Corridor)+
  ggtitle("Shannon diversity last day")+
  theme_bw()

###all together
ggplot(data_lastday, aes( x = Temp_Regime, y = Sh_div, col = Temp_Regime))+
  geom_boxplot()+
  ggtitle("Shannon diversity last day")+
  theme_bw()



##we use the normal error dist to begin with identity link
#code the model
mod_last_day_Sh_1_ID<-glm(Sh_div~ Corridor*Temp_Regime,
                          data = data_lastday, family = gaussian(link = "identity"))
##visual checking of the model
plot(mod_last_day_Sh_1_ID)

##look at the R squared to have an idea of how much variance is explained 
rsquared(mod_last_day_Sh_1_ID)

###plot the residual of the model to look at their distribution

hist(resid(mod_last_day_Sh_1_ID))

###we test the normality of the residuals
shapiro.test(resid(mod_last_day_Sh_1_ID))
##it is good! Looks like they are normal

###dharma testing is also happy with it
sim5<-simulateResiduals(mod_last_day_Sh_1_ID, n = 1000)
plot(sim5)


##model checking

anova(mod_last_day_Sh_1_ID,test = "F")
Anova(mod_last_day_Sh_1_ID)
summary(mod_last_day_Sh_1_ID)

##post hoc test
summary(glht(mod_last_day_Sh_1_ID, mcp(Temp_Regime = "Tukey")))


####doin the simple ANOVA without model
plot(lm(Sh_div~ Corridor*Temp_Regime,
       data = data_lastday))

data_lastday %>% 
  group_by(Temp_Regime, Corridor) %>%
  shapiro_test(Sh_div)

ggqqplot(data_lastday, "Sh_div", ggtheme = theme_bw()) +
  facet_grid(Temp_Regime ~ Corridor)

data_lastday %>% levene_test(Sh_div~ Corridor*Temp_Regime)

anova_test(Sh_div~ Corridor*Temp_Regime,
           data = data_lastday)

data_lastday %>% pairwise_t_test(
  Sh_div ~ Temp_Regime, 
  p.adjust.method = "bonferroni"
)

dd<-data_lastday %>% 
  group_by(Temp_Regime) %>% 
  dplyr::summarise(sd = sd(Sh_div))



#####BETA DIVERSITY BETWEEN THE TWO PATCHES ACROSS CORRIDOR AND TREATMENTS at the end of the experiment####

#####BETA DIVERSITY CALCULATION#################
data_Beta<-dataG%>%
  dplyr::select(Date, Replicate, Temp_Regime, Corridor,Patch, Total_Blepharisma,
                Total_Didinium, Total_Colpidium, Total_Homolazoon, Total_Bursaria,
                Total_Spirostomum, Total_Paramecium) %>% 
  dplyr::ungroup() %>% 
  group_by(Date,Replicate,Temp_Regime,Corridor)


###calculate the Bdiv indeces for each replica of each treatment

data_Beta <- data_Beta %>% group_by(Date,Replicate,Temp_Regime,Corridor) %>%
  tidyr::nest() %>%
  mutate(data = purrr::map(data, ~.x %>% dplyr::select(-Patch))) %>%
  mutate(data = purrr::map(data, ~.x %>%
                             mutate(beta = betapart::beta.multi.abund(.x)$beta.BRAY,
                                    Patch = c("A","B"))))%>%
  tidyr::unnest(cols = c(data))


data_Beta<-data_Beta %>% dplyr::ungroup() %>%
  mutate(NumDays = as.numeric(difftime(data_Beta$Date, data_Beta$Date[1], units = "days")))

ggplot(data_Beta, aes(x=Date, y =beta,
                            col= Temp_Regime, group = Replicate , fill=Temp_Regime ))+
  geom_point()+
  geom_smooth(aes(group=Temp_Regime),alpha = 0.1)+
  facet_wrap(~Corridor)+
  ggtitle("Patch-dissimilarity index (Beta Bray)")+
  theme_bw()

ggplot(data_Beta, aes(x=Date, y =beta, color = Corridor))+
  geom_point()+
  geom_smooth(alpha = 0.1)+
  ggtitle("Patch-dissimilarity index (Beta Bray)")+
  theme_bw()

###select just the last day to analyise the Beta div lst day
dataBeta_lastday<-data_Beta %>% filter(Date =="2021-12-20") %>% filter((Patch=="A"))

#we need factors to avoid problems in the model
dataBeta_lastday$Temp_Regime<-as.factor(dataBeta_lastday$Temp_Regime)
dataBeta_lastday$Corridor<-as.factor(dataBeta_lastday$Corridor)
write.csv(dataBeta_lastday, file="Data_Betadiv_lastday.csv")

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

ggplot(dataBeta_lastday, aes(x = Corridor, y = beta, col = Corridor))+
  geom_boxplot()+
  ggtitle("Beta diversity between two patches LAST DAY")+
  theme_bw()

##plot the data to see distribution
ggplot(dataBeta_lastday, aes(x = beta))+
  geom_histogram(bins = 30,position="identity", alpha = .5)+
  facet_grid(~Corridor)+
  ggtitle("Beta div between patches last Day")+
  theme_bw()

ggplot(dataBeta_lastday, aes(x = beta))+
  geom_histogram(bins = 30,position="identity", alpha = .5)+
  ggtitle("Beta div between patches last Day")+
  theme_bw()


#model WITH interaction####

mod_last_day_Beta_INT<-glm(beta~ Corridor*Temp_Regime,
                           data = dataBeta_lastday, family = "gaussian")

mod_last_day_Beta_gamma<-glm(beta+1~ Corridor*Temp_Regime,
                            data = dataBeta_lastday, family = Gamma(link = "identity"))

##visual checking of the model
plot(mod_last_day_Beta_INT)
plot(mod_last_day_Beta_gamma)


simulateResiduals(mod_last_day_Beta_INT, n = 1000, plot = T)
simulateResiduals(mod_last_day_Beta_gamma, n = 1000, plot = T)

plot(res4)
##look at the R squared to have an idea of how much variance is explained 
rsquared(mod_last_day_Beta_INT)
rsquared(mod_last_day_Beta_gamma)

###plot the residual of the model to look at their distribution
hist(resid(mod_last_day_Beta_INT))

hist(resid(mod_last_day_Beta_gamma))

###we test the normality of the residuals
shapiro.test(resid(mod_last_day_Beta_INT))
shapiro.test(resid(mod_last_day_Beta_gamma))


summary(mod_last_day_Beta_INT)
summary(mod_last_day_Beta_gamma)


Anova(mod_last_day_Beta_gamma)
##post hoc test
summary(glht(mod_last_day_Beta_INT_SQ, mcp(Temp_Regime = "Tukey")))

#####

####Analysis with the simple two way ANOVA without model#####
dataBeta_lastday %>% ungroup()


modB<-lm(sqrt(beta)~ Corridor*Temp_Regime,
         data = dataBeta_lastday)

plot(modB)
hist(residuals(modB))


shapiro.test(residuals(modB))

dataBeta_lastday<-dataBeta_lastday %>% 
  group_by(Temp_Regime, Corridor) %>%
  dplyr::mutate(sq_beta = sqrt(beta)) 
  
dataBeta_lastday %>% shapiro_test(sq_beta)

ggqqplot(dataBeta_lastday, "sq_beta", ggtheme = theme_bw()) +
  facet_grid(Temp_Regime ~ Corridor)

dataBeta_lastday %>% ungroup() %>%  
  levene_test(sq_beta ~ Corridor*Temp_Regime)

dataBeta_lastday<-dataBeta_lastday %>% ungroup()

anova_test(sq_beta ~ Corridor*Temp_Regime,
           data = dataBeta_lastday)

summary(aov(modB))


pairwise_t_test(dataBeta_lastday, sq_beta ~ Temp_Regime, 
  p.adjust.method = "bonferroni"
)








#############Analysis THROUGH TIME####################


#####ANALYSIS of SPECIES RICHNESS THROUHG TIME--------------
###create a columns with the interaction factor in case we need it for post hoc comparisons

##plot the data to see distribution
ggplot(data_Sh_div, aes(x = Sp_rich))+
  geom_histogram(bins = 30,position="identity", alpha = .5)+
  ggtitle("Species richness through time")+
  theme_bw()

#with poisson model ~ corridors

ggplot(data_Sh_div %>% group_by(Temp_Regime, Corridor,Replicate),
       aes(x = NumDays, y = Sp_rich,  fill = Temp_Regime))+
  geom_point(shape = 21, alpha = 0.3, size = 2,
             position=position_jitter(width=0.3, height=0.05, seed = 1))+
  geom_smooth(method = "loess",aes(color = Temp_Regime),se = T, alpha = 0.2)+
  facet_grid(~Corridor)+
  ggtitle("Species Richness")+
  theme_bw()


ggplot(data_Sh_div %>% group_by(Temp_Regime, Corridor,Replicate),
       aes(x = Temp_Regime, y = Sp_rich,  col = Temp_Regime))+
  geom_boxplot()+
  theme_bw()
  



#run the function to detect which distribution the data fits

descdist(data_Sh_div$Sp_rich, discrete = FALSE, boot = 500)

colnames(data_Sh_div)

###code the model

#model with all interactions with TIME
##create a variable with the initial richness to use in the model to fit the intercept
data_Sh_div$initial_richness<-data_Sh_div$Sp_rich[1]

mod_Sp_time_NOINT<-glm(Sp_rich~NumDays +
                   Temp_Regime + Corridor +NumDays:Temp_Regime + NumDays:Corridor + Corridor:Temp_Regime+
                   Temp_Regime:Corridor:NumDays,
                 data = data_Sh_div,
                 family = "poisson" )


##visual checking of the model
plot(mod_Sp_time_NOINT)
plot(mod_simp)

##look at the R squared to have an idea of how much variance is explained 
rsquared(mod_Sp_time_NOINT)

###plot the residual of the model to look at their distribution

hist(resid(mod_Sp_time_NOINT))
shapiro.test(resid(mod_Sp_time_NOINT))


##model checking
summary(mod_Sp_time_NOINT)
anova(mod_Sp_time_NOINT,test = "Chisq")
Anova(mod_Sp_time_NOINT)
Anova(mod_Sp_time_NOINT, test.statistic = "F")
anova(mod_Sp_time_NOINT,test = "F")
##post hoc test
summary(glht(mod_Sp_time_NOINT, mcp(Temp_Regime = "Tukey")))


####model simplification
mod_simp1<-glm(Sp_rich~NumDays +
                 Temp_Regime + Corridor +NumDays:Temp_Regime + NumDays:Corridor + Corridor:Temp_Regime,
               data = data_Sh_div,
              family = "poisson" )


##visual checking of the model
plot(mod_simp1)

##look at the R squared to have an idea of how much variance is explained 
rsquared(mod_simp1)

###plot the residual of the model to look at their distribution

hist(resid(mod_simp1))
shapiro.test(resid(mod_simp1))


##model checking
summary(mod_simp1)
anova(mod_simp1,test = "Chisq")
Anova(mod_simp1)


####simplest model##
mod_simp2<-glm(Sp_rich~NumDays +
                 Temp_Regime + Corridor,
               data = data_Sh_div,
               family = "poisson" )


##visual checking of the model
plot(mod_simp2)

##look at the R squared to have an idea of how much variance is explained 
rsquared(mod_simp2)

###plot the residual of the model to look at their distribution

hist(resid(mod_simp2))
shapiro.test(resid(mod_simp2))


##model checking
summary(mod_simp2)
anova(mod_simp2,test = "Chisq")
Anova(mod_simp2)

##post hoc
summary(glht(mod_simp2, mcp(Temp_Regime = "Tukey")))

####simplest model PART 2##
mod_simp3<-glm(Sp_rich~NumDays +
                 Temp_Regime,
               data = data_Sh_div,
               family = "poisson" )

mod_simp4<-glm(Sp_rich~NumDays,
               data = data_Sh_div,
               family = "poisson" )

##visual checking of the modell
plot(mod_simp3)

##look at the R squared to have an idea of how much variance is explained 
rsquared(mod_simp3)

###plot the residual of the model to look at their distribution

hist(resid(mod_simp3))
shapiro.test(resid(mod_simp3))


##model checking
summary(mod_simp3)
anova(mod_simp3,test = "Chisq")
rstatix::Anova(mod_simp3)
Anova(mod_simp3)
Anova(mod_simp3, type = 3)
##post hoc
summary(glht(mod_simp3, mcp(Temp_Regime = "Tukey")))


lrtest(mod_Sp_time_NOINT,mod_simp1)
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
  ggtitle("Poisson model without intercept")+
  theme_bw()



#model SPECIES RICHNESS TIME using splines####

##create a variable with the initial richness to use in the model to fit the intercept
data_Sh_div$initial_richness<-data_Sh_div$Sp_rich[1]

mod_Sp_time_SPLINES<-glm(Sp_rich~ns(NumDays,3) +
                   Temp_Regime + Corridor +ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor +
                   Temp_Regime:Corridor:ns(NumDays,3),
                 data = data_Sh_div,
                 family = "poisson" )

##visual checking of the model
plot(mod_Sp_time_SPLINES)


###plot the residual of the model to look at their distribution

hist(resid(mod_Sp_time_SPLINES))
shapiro.test(resid(mod_Sp_time_SPLINES))

simulateResiduals(mod_Sp_time_SPLINES, plot = T)

#model checking
summary(mod_Sp_time_SPLINES)
anova(mod_Sp_time_SPLINES,test = "Chisq")
Anova(mod_Sp_time_SPLINES)


##post hoc test
summary(glht(mod_Sp_time_SPLINES, mcp(Temp_Regime = "Tukey")))


##predict values from the model and plot it, to see the fit of the model to real data
gg.df_POI_SPLI <- data.frame(preds = predict(mod_Sp_time_SPLINES, newdata = data_Sh_div,type = "response" ),
                        upr = predict(mod_Sp_time_SPLINES, newdata = data_Sh_div ,se.fit = T,type = "response" )$fit +
                          predict(mod_Sp_time_SPLINES, newdata = data_Sh_div,se.fit = T,type = "response" )$se.fit ,
                        lwr = predict(mod_Sp_time_SPLINES, newdata = data_Sh_div ,se.fit = T ,type = "response")$fit -
                          predict(mod_Sp_time_SPLINES, newdata = data_Sh_div ,se.fit = T ,type = "response")$se.fit , data_Sh_div)


ggplot(gg.df_POI_SPLI, aes(x =NumDays, y= Sp_rich,
                      col =Temp_Regime, fill = Temp_Regime))+
  geom_point(position=position_jitter(width=0.3, height=0.05, seed = 1), alpha = .2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),col="transparent",alpha=0.1)+
  facet_wrap(~Corridor)+
  ggtitle("Poisson model with intercept")+
  theme_bw()

###model with splines is visually better but it does not fit properly (residuals are not significantly normal according to shapiro)






##

#####ANALYSIS of SHANNON DIVERSITY THROUHG TIME#########

##plot the data to see distribution
data_Sh_div %>% group_by(Replicate, Corridor, Temp_Regime)
ggplot(data_Sh_div, aes(x = Sh_div))+
  geom_histogram(bins = 30,position="identity", alpha = .5)+
  ggtitle("Shannon Diversity through time")+
  theme_bw()


ggplot(data_Sh_div, aes(x = NumDays, y = Sh_div, fill = Temp_Regime),col = Temp_Regime,
       group = Replicate)+
  geom_point(aes(col=Temp_Regime), alpha=0.2)+
  geom_smooth(aes(col=Temp_Regime),alpha=0.2)+
  ggtitle("Shannon Diversity through time")+
  facet_grid(~Corridor)+
  theme_bw()


###code the model with cubic Spline to follow the trend through time
###create the variable with the initial value of Shannon diversity at day 0 to be used as baseline
data_Sh_div$initial_Sh<-data_Sh_div$Sh_div[1]

mod_SH_time1_SPLI_noint<-glm(Sh_div~ ns(NumDays,3) +
                         Temp_Regime + 
                         Corridor + 
                         ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor + Temp_Regime:Corridor+
                         ns(NumDays,3):Temp_Regime:Corridor, data =data_Sh_div, family = "gaussian")

##visual checking of the model
plot(mod_SH_time1_SPLI_noint)   


hist(residuals(mod_SH_time1_SPLI_noint))
shapiro.test(residuals(mod_SH_time1_SPLI_noint))

###results of the model

summary(mod_SH_time1_SPLI_noint)
Anova(mod_SH_time1_SPLI_noint, test.statistic = "F")
Anova(mod_SH_time1_SPLI_noint, type = 3)
anova(mod_SH_time1_SPLI_noint, test = "F")


####no interaction of the cubic term, we take it out of the model and re run

mod_SH_time1_SPLI_Simp1<-glm(Sh_div~ ns(NumDays,3) +
                               Temp_Regime + 
                               Corridor + 
                               ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor + Temp_Regime:Corridor,
                             data =data_Sh_div, family = "gaussian")

##visual checking of the model
plot(mod_SH_time1_SPLI_Simp1)   

hist(residuals(mod_SH_time1_SPLI_Simp1))
shapiro.test(residuals(mod_SH_time1_SPLI_Simp1))

###results of the model

summary(mod_SH_time1_SPLI_Simp1)
Anova(mod_SH_time1_SPLI_Simp1)
anova(mod_SH_time1_SPLI_Simp1, test = "F")
Anova(mod_SH_time1_SPLI_Simp1, test.statistic = "F")


####no effect of the Temp_Regime:Corridor interaction, we take it out of the model and re run

mod_SH_time1_SPLI_Simp2<-glm(Sh_div~ ns(NumDays,3) +
                               Temp_Regime + 
                               Corridor + 
                               ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor,
                             data =data_Sh_div, family = "gaussian")

unique(data_Sh_div$Sh_div)
##visual checking of the model
plot(mod_SH_time1_SPLI_Simp2)   

hist(residuals(mod_SH_time1_SPLI_Simp2))
shapiro.test(residuals(mod_SH_time1_SPLI_Simp2))

summary(glht(mod_SH_time1_SPLI_Simp2, mcp(Temp_Regime = "Tukey")))


###results of the model

summary(mod_SH_time1_SPLI_Simp2)
Anova(mod_SH_time1_SPLI_Simp2, test.statistic = "F", type = 3)
Anova(mod_SH_time1_SPLI_Simp2, test.statistic = "F", type = 2)
anova(mod_SH_time1_SPLI_Simp2, test = "F")


####likelyhood ratio test to compare the model and justify the step wise model selection
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


###fit a model just with the significant interaction of corridor with time to 
##display the differences
modCor<-glm(Sh_div~ ns(NumDays,3)+Corridor + 
  + ns(NumDays,3):Corridor,
data =data_Sh_div, family = "gaussian")

plot(modCor)
hist(residuals(modCor))
shapiro_test(residuals(modCor))

summary(modCor)

gg.df_CORR <- data.frame(preds = predict(modCor, newdata = data_Sh_div,type = "response" ),
                          upr = predict(modCor, newdata = data_Sh_div ,se.fit = T,type = "response" )$fit +
                            predict(modCor, newdata = data_Sh_div,se.fit = T,type = "response" )$se.fit ,
                          lwr = predict(modCor, newdata = data_Sh_div ,se.fit = T ,type = "response")$fit -
                            predict(modCor, newdata = data_Sh_div ,se.fit = T ,type = "response")$se.fit , data_Sh_div)


ggplot(gg.df_CORR, aes(x =NumDays, y= Sh_div,group = Corridor,
                       col =Corridor, fill = Corridor))+
  geom_point(alpha = 0.2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  theme_bw()

###fit a model just with the significant interaction of TEMP TEGIME with time to 
##display the differences
modTemp<-glm(Sh_div~ ns(NumDays,3)+Temp_Regime + 
              + ns(NumDays,3):Temp_Regime,
            data =data_Sh_div, family = "gaussian")

plot(modTemp)
hist(residuals(modTemp))
shapiro_test(residuals(modTemp))

summary(modTemp)

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







####ANALYSIS of Shannon diversity VARIANCE across replicate through time#####
#calculating median Sh_diversity between replicates + upper and lower standard deviations + variance

data_Sh_div_VAR<-data_Sh_div %>%
  dplyr::group_by(Date,Temp_Regime, Corridor, NumDays) %>% 
  dplyr::summarise(median.Sh_div = median(Sh_div), sd.val = sd(Sh_div), Sh_variance = var(Sh_div),
                   upr = median.Sh_div + sd.val, lwr = median.Sh_div - sd.val)


ggplot(data_Sh_div_VAR, aes(x = Date, y = Sh_variance   , color = Temp_Regime,
                            fill=Temp_Regime,))+
  geom_point(size = 2, alpha = 0.4)+
  geom_smooth(method = "loess",alpha= 0.2)+
  facet_grid(~Corridor)+
  ggtitle("Shannon diversity variance")+
  theme_bw()


####avoid condiering the first day in which we have variance zero
data_Sh_div_VAR_no_ZERO<-data_Sh_div_VAR %>%filter(NumDays != 0)

data_Sh_div_VAR_no_ZERO$Temp_Regime<-as.factor(data_Sh_div_VAR_no_ZERO$Temp_Regime)
data_Sh_div_VAR_no_ZERO$Corridor<-as.factor(data_Sh_div_VAR_no_ZERO$Corridor)

str(data_Sh_div_VAR_no_ZERO)
ggplot(data_Sh_div_VAR_no_ZERO, aes(x = Date, y = Sh_variance   , color = Temp_Regime,
                            fill=Temp_Regime,))+
  geom_point(size = 2, alpha = 0.4)+
  geom_smooth(method = "loess",alpha= 0.2)+
  facet_grid(~Corridor)+
  ggtitle("Shannon diversity variance")+
  theme_bw()

##plot the data to see distribution
ggplot(data_Sh_div_VAR_no_ZERO, aes(x = Sh_variance))+
  geom_histogram(bins = 30,position="identity", alpha = .5)+
  ggtitle("Shannon Diversity variance through time")+
  theme_bw()


##code the model with splines 1
mod_SH_VAR_SPL1<-glm(Sh_variance~ ns(NumDays,3) +
                               Temp_Regime + 
                               Corridor + 
                               ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor + Temp_Regime:Corridor+
                               ns(NumDays,3):Temp_Regime:Corridor, data =data_Sh_div_VAR_no_ZERO,
                     family = gaussian(link = "identity"))



##visual checking of the model
plot(mod_SH_VAR_SPL1)   
hist(residuals(mod_SH_VAR_SPL1))
shapiro.test(residuals(mod_SH_VAR_SPL1))

###results of the model

summary(mod_SH_VAR_SPL1)
Anova(mod_SH_VAR_SPL1, test.statistic = "F")
Anova(mod_SH_VAR_SPL1, type = 3)
anova(mod_SH_VAR_SPL1, test = "F")



##post hoc
anova_var<-aov(Sh_variance~
      Temp_Regime*Corridor,
      data =data_Sh_div_VAR_no_ZERO)
summary(anova_var)

sig_comp<-as.data.frame(?TukeyHSD(anova_var)$"Temp_Regime:Corridor")
sig_comp<-sig_comp %>% filter(`p adj` < 0.05)


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
  ggtitle("Model with cubic interaction")+
  theme_bw()

ggplot(data_Sh_div_VAR_no_ZERO, aes(x = Temp_Regime, y = Sh_variance, col = Corridor))+
  geom_boxplot()+
  theme_bw()

ggplot(data_Sh_div_VAR_no_ZERO, aes(x = Corridor, y = Sh_variance, col = Corridor))+
  geom_boxplot()+
  theme_bw()


##code the simpler model with splines####
mod_SH_VAR_SPL2<-glm(Sh_variance~ ns(NumDays,3) +
                       Temp_Regime + 
                       Corridor + 
                       ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor + Temp_Regime:Corridor,
                     data =data_Sh_div_VAR_no_ZERO,
                     family = gaussian(link = "identity"))

##visual checking of the model
plot(mod_SH_VAR_SPL2)   
hist(residuals(mod_SH_VAR_SPL2))
shapiro.test(residuals(mod_SH_VAR_SPL2))

###results of the model

summary(mod_SH_VAR_SPL2)
Anova(mod_SH_VAR_SPL2, test.statistic = "F")
Anova(mod_SH_VAR_SPL2, type = 3)
anova(mod_SH_VAR_SPL2, test = "F")

##predict values from the model and plot it, to see the fit of the model to real data
gg.df_VAR2 <- data.frame(preds = predict(mod_SH_VAR_SPL2, newdata = data_Sh_div_VAR_no_ZERO,type = "response" ),
                        upr = predict(mod_SH_VAR_SPL2, newdata = data_Sh_div_VAR_no_ZERO ,se.fit = T,type = "response" )$fit +
                          predict(mod_SH_VAR_SPL2, newdata = data_Sh_div_VAR_no_ZERO,se.fit = T,type = "response" )$se.fit ,
                        lwr = predict(mod_SH_VAR_SPL2, newdata = data_Sh_div_VAR_no_ZERO ,se.fit = T ,type = "response")$fit -
                          predict(mod_SH_VAR_SPL2, newdata = data_Sh_div_VAR_no_ZERO ,se.fit = T ,type = "response")$se.fit , data_Sh_div_VAR_no_ZERO)


ggplot(gg.df_VAR2, aes(x =NumDays, y= Sh_variance,group = Temp_Regime, col =Temp_Regime, fill = Temp_Regime))+
  geom_point(alpha = 0.2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  facet_grid(~Corridor)+
  ggtitle("Model WITHOUT cubic interaction")+
  theme_bw()


####comparing the two models
AIC(mod_SH_VAR_SPL1, mod_SH_VAR_SPL2)
lrtest(mod_SH_VAR_SPL2,mod_SH_VAR_SPL1)
###

####
######ANALYSIS of BETA DIVERSITY between the two patches THROUHG TIME#########

data_Beta<-data_Beta %>% group_by(Date, Replicate, Corridor, Temp_Regime, Patch) %>%
  filter(Patch == "A")


###plot the data with loess to see the trends through time
ggplot(data_Beta,aes(x = NumDays, y = beta, col = Temp_Regime, fill = Temp_Regime))+
  geom_point(aes(col =Temp_Regime), alpha = 0.2)+
  geom_smooth(method = "loess", alpha =  0.2)+
  facet_grid(~Corridor)+
  theme_bw()

##histogram of the data to see distribution
data_Beta %>% group_by(Replicate, Corridor, Temp_Regime)

ggplot(data_Beta, aes(x = beta))+
  geom_histogram(bins = 30,position="identity", alpha = .5)+
  facet_grid(~Corridor)+
  ggtitle("Patch beta Diversity through time")+
  theme_bw()


ggplot(data_Beta_no_zero, aes(x = beta))+
  geom_histogram(bins = 30,position="identity", alpha = .5)+
  ggtitle("Patch beta Diversity through time")+
  theme_bw()

###model with splines
library(betareg)

###take away the first day from the dataset
data_Beta_no_zero<-data_Beta %>% filter(NumDays != 0)
ggplot(data_Beta_no_zero,aes(x = NumDays, y = beta, col = Temp_Regime, fill = Temp_Regime))+
  geom_point(aes(col =Temp_Regime), alpha = 0.2)+
  geom_smooth(method = "loess", alpha =  0.2)+
  facet_grid(~Corridor)+
  theme_bw()

mod_beta_time_SPLI<-glm(beta ~ ns(NumDays,3) +
                               Temp_Regime + 
                               Corridor + 
                               ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor +
                               ns(NumDays,3):Temp_Regime:Corridor, data =data_Beta_no_zero, 
                        family = gaussian(link = "identity"))



data_Beta_no_zero<-data_Beta_no_zero %>% mutate(beta_trans = ((beta*(702-1)+0.5))/702)


mod_beta_time_SPLI_beta<-betareg(beta_trans ~ ns(NumDays,3) +
                                   Temp_Regime + 
                                   Corridor + 
                                   ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor +
                                   ns(NumDays,3):Temp_Regime:Corridor, data =data_Beta_no_zero)
mod_beta_time_SPLI_beta<-betareg(beta_trans ~ ns(NumDays,3) +
                                   Temp_Regime + 
                                   Corridor + 
                                   ns(NumDays,3):Temp_Regime + ns(NumDays,3):Corridor, data =data_Beta_no_zero)
##visual checking of the model
plot(mod_beta_time_SPLI)   

hist(residuals(mod_beta_time_SPLI))
shapiro.test(residuals(mod_beta_time_SPLI))
shapiro_test(residuals(mod_beta_time_SPLI))



plot(mod_beta_time_SPLI_beta)   
plot(mod_beta_time_SPLI_beta, which = 1:4, type = "pearson")
plot(mod_beta_time_SPLI_beta, which = 5, type = "deviance", sub.caption = "")
plot(mod_beta_time_SPLI_beta, which = 1, type = "deviance", sub.caption = "")
##dharma checking
simulateResiduals(mod_beta_time_SPLI, n = 1000, plot = T)


###results of the model
summary(mod_beta_time_SPLI)     

Anova(mod_beta_time_SPLI)
Anova(mod_beta_time_SPLI, type = 3)

summary(mod_beta_time_SPLI_beta)     

Anova(mod_beta_time_SPLI_beta)

##predict values from the model and plot it, to see the fit of the model to real data
gg.df.beta <- data.frame(preds = predict(mod_beta_time_SPLI, newdata = data_Beta,type = "response" ),
                    upr = predict(mod_beta_time_SPLI, newdata = data_Beta ,se.fit = T,type = "response" )$fit +
                      predict(mod_beta_time_SPLI, newdata = data_Beta,se.fit = T,type = "response" )$se.fit ,
                    lwr = predict(mod_beta_time_SPLI, newdata s= data_Beta ,se.fit = T ,type = "response")$fit -
                      predict(mod_beta_time_SPLI, newdata = data_Beta ,se.fit = T ,type = "response")$se.fit , data_Beta)


ggplot(gg.df.beta, aes(x =NumDays, y= beta+1,col =Temp_Regime, fill = Temp_Regime))+
  geom_point()+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  facet_grid(~Corridor)+
  theme_bw()



gg.df.beta.beta <- data.frame(preds = betareg::predict(mod_beta_time_SPLI_beta, newdata = data_Beta,type = "response" ),
                              upp= betareg::predict(mod_beta_time_SPLI_beta, newdata = data_Beta,type = "quantile", at = 0.95),
                              low=betareg::predict(mod_beta_time_SPLI_beta, newdata = data_Beta,type = "quantile", at = 0.05), data_Beta)


ggplot(gg.df.beta.beta, aes(x =NumDays, y= beta_trans,col =Temp_Regime,
                            fill = Temp_Regime))+
  geom_point()+
  geom_line(aes(y=preds))+
#geom_ribbon(aes(ymin = low,
               # ymax = upp),colour="transparent", alpha=0.2)+
  facet_grid(~Corridor)+
  theme_bw()


gg.df.beta.NOSPL <- data.frame(preds = betareg::predict(mod_beta_time_beta, newdata = data_Beta,type = "response" ),
                              upp= betareg::predict(mod_beta_time_beta, newdata = data_Beta,type = "quantile", at = 0.95),
                              low=betareg::predict(mod_beta_time_beta, newdata = data_Beta,type = "quantile", at = 0.05), data_Beta)


ggplot(gg.df.beta.NOSPL, aes(x =NumDays, y= beta_trans,col =Temp_Regime,
                            fill = Temp_Regime))+
  geom_point()+
  geom_line(aes(y=preds))+
  #geom_ribbon(aes(ymin = low,
  # ymax = upp),colour="transparent", alpha=0.2)+
  facet_grid(~Corridor)+
  theme_bw()
#####end####