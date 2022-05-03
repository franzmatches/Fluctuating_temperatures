#########################################################################################
### ANALYSIS THROUGH TIME ###
#########################################################################################

### Preamble ###
library(ggh4x)
library(tidyverse) # data wrangling and plotting
library(ggpubr) # multi panel figures
library(rstatix) # pipe-friendly model assumption tests
library(multcomp) # parametric model comparison functions
library(lmtest) # model comparisons

source("Code/data_preparation.R")


#CODE FOR PLOTTING INDIVIDUAL SPECIES THROUGH TIME AT LANSCAPE AND PATCH LEVEL

#pivot for obtaining the species variable for plotting
data_pivot<-data_patches%>%
  pivot_longer(cols = -c(Date, Replicate, Temp_Regime,
                         Corridor,Patch, beta, NumDays),
               names_to ="species",
               values_to= "abundance")%>%
  group_by(Temp_Regime, Corridor,Patch,Replicate)


data_merged_corridors<-data_pivot %>% ungroup () %>% 
  group_by(Date, Replicate, Temp_Regime,Patch, NumDays, species) %>% 
  summarize(summed_abundance = sum(abundance))


ggplot(data_pivot, aes(x = NumDays, y = log(abundance+1), 
                       colour = species, fill = species))+
  geom_line(aes(linetype = (as.factor(Replicate))), size = 0.7)+
  facet_nested(Corridor + Patch ~ Temp_Regime,scales = "fixed",
               labeller = label_value,
               strip = ggh4x::strip_nested(size="constant",bleed=T))+
  theme_bw()

ggsave("Log_species abundance patch and corridor_replicates.pdf", units="in", width=16, height=10)


ggplot(data_pivot, aes(x = NumDays, y = log(abundance+1), 
                       colour = species, fill = species))+
  geom_smooth(method = "loess",alpha = 0.15)+
  facet_nested(Corridor + Patch ~ Temp_Regime,scales = "fixed",
               labeller = label_value,
               strip = ggh4x::strip_nested(size="constant",bleed=T))+
  theme_bw()

ggsave("Log_species abundance patch and corridor_smooth.pdf", units="in", width=16, height=10)


#----------------------------------------------------------------------------------------
## ANALYSIS of SPECIES RICHNESS THROUGH TIME ##
#----------------------------------------------------------------------------------------

##model with all interactions with the time factor (in this case NumDays:number of days of the experiment)
mod_Sp_time<-glm(Sp_rich~NumDays +
                   Temp_Regime + Corridor +NumDays:Temp_Regime + NumDays:Corridor + Corridor:Temp_Regime+
                   Temp_Regime:Corridor:NumDays,
                 data = data_pooled,
                 family = "poisson" )


##visual checking of the model
plot(mod_Sp_time)

##plot the residual of the model to look at their distribution
hist(resid(mod_Sp_time))

##check normality of the residuals with Shapiro-Wilk test
rstatix::shapiro_test(resid(mod_Sp_time))

##model outputs
summary(mod_Sp_time)
Anova(mod_Sp_time)


##model simplification 1 (without three way interaction)
mod_simp1<-glm(Sp_rich~NumDays +
                 Temp_Regime + Corridor +NumDays:Temp_Regime + NumDays:Corridor + Corridor:Temp_Regime,
               data = data_pooled,
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


##model simplification 2 (without interactions)
mod_simp2<-glm(Sp_rich~NumDays +
                 Temp_Regime + Corridor,
               data = data_pooled,
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
               data = data_pooled,
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
summary(multcomp::glht(mod_simp3, multcomp::mcp(Temp_Regime = "Tukey")))

##control model just with time to be sure of the significant effect of the Temperature Regime Factor
mod_simp4<-glm(Sp_rich~NumDays,
               data = data_pooled,
               family = "poisson" )

##Likelihood Ratio Test to compare if our model selection is significant 
lmtest::lrtest(mod_Sp_time,mod_simp1)
lmtest::lrtest(mod_simp1,mod_simp2)
lmtest::lrtest(mod_simp2,mod_simp3)
lmtest::lrtest(mod_simp3,mod_simp4)

##predict values from the model and plot it, to see the fit of the model to real data
df.simp<-data.frame(preds = predict(mod_simp3, newdata = data_pooled,type = "response" ),
                    upr = predict(mod_simp3, newdata = data_pooled ,se.fit = T,type = "response" )$fit + predict(mod_simp3, newdata = data_pooled,se.fit = T,type = "response" )$se.fit ,
                    lwr = predict(mod_simp3, newdata = data_pooled ,se.fit = T ,type = "response")$fit - predict(mod_simp3, newdata = data_pooled ,se.fit = T ,type = "response")$se.fit , 
                    data_pooled)

ggplot(df.simp, aes(x =NumDays, y= Sp_rich,
                    col =Temp_Regime, fill = Temp_Regime))+
  geom_point(position=position_jitter(width=0.3, height=0.05, seed = 1), alpha = .2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),col="transparent",alpha=0.1)+
  #facet_wrap(~Corridor)+
  ggtitle("Species Richness Poisson model")+
  theme_bw()

#----------------------------------------------------------------------------------------
## ANALYSIS of SHANNON DIVERSITY THROUHG TIME ##
#----------------------------------------------------------------------------------------

####Build the the model with cubic Spline for the time variable 
##model with all interactions with the time factor (in this case NumDays:number of days of the experiment)
mod_SH_time1_SPLI<-glm(Sh_div~ splines::ns(NumDays,3) +
                         Temp_Regime + 
                         Corridor + 
                         splines::ns(NumDays,3):Temp_Regime + splines::ns(NumDays,3):Corridor + Temp_Regime:Corridor+
                         splines::ns(NumDays,3):Temp_Regime:Corridor, data =data_pooled, family = "gaussian")

##visual checking of the model
plot(mod_SH_time1_SPLI)   

##plot the residual of the model to look at their distribution
hist(residuals(mod_SH_time1_SPLI))

##check normality of the residuals with Shapiro-Wilk test
shapiro.test(residuals(mod_SH_time1_SPLI))

##Output of the model
summary(mod_SH_time1_SPLI)
Anova(mod_SH_time1_SPLI, test.statistic = "F")

##model simplification 1 - no three way interaction
mod_SH_time1_SPLI_Simp1<-glm(Sh_div~ splines::ns(NumDays,3) +
                               Temp_Regime + 
                               Corridor + 
                               splines::ns(NumDays,3):Temp_Regime + splines::ns(NumDays,3):Corridor + Temp_Regime:Corridor,
                             data =data_pooled, family = "gaussian")

##visual checking of the model
plot(mod_SH_time1_SPLI_Simp1)   

##plot the residual of the model to look at their distribution
hist(residuals(mod_SH_time1_SPLI_Simp1))

##check normality of the residuals with Shapiro-Wilk test
shapiro.test(residuals(mod_SH_time1_SPLI_Simp1))

###output of the model
summary(mod_SH_time1_SPLI_Simp1)
Anova(mod_SH_time1_SPLI_Simp1, test.statistic = "F")

##model simplification 2 
#no effect of the Temp_Regime:Corridor interaction, we take it out of the model and re run
mod_SH_time1_SPLI_Simp2<-glm(Sh_div~ splines::ns(NumDays,3) +
                               Temp_Regime + 
                               Corridor + 
                               splines::ns(NumDays,3):Temp_Regime + splines::ns(NumDays,3):Corridor,
                             data =data_pooled, family = "gaussian")


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

ggplot(data_pooled, aes(x = Temp_Regime, y = Sh_div, col = Temp_Regime))+
  geom_boxplot()+
  ggtitle("Shannon diversity all time points together")+
  theme_bw()

###visually it does not look like we have a large difference between temperature regimes 
# on average 
##running an anova to check  that there is no significant difference
aov_mod_SH_time1_SPLI_Simp2<-aov(Sh_div~
                                   Temp_Regime,
                                 data =data_pooled)
summary(aov_mod_SH_time1_SPLI_Simp2)


##Likelihood Ratio Test to compare if our model selection is significant 
lmtest::lrtest(mod_SH_time1_SPLI,mod_SH_time1_SPLI_Simp1)
lmtest::lrtest(mod_SH_time1_SPLI_Simp1,mod_SH_time1_SPLI_Simp2)

##predict values from the model and plot it, to see the fit of the model to real data
gg.df_noint <- data.frame(preds = predict(mod_SH_time1_SPLI_Simp2, newdata = data_pooled,type = "response" ),
                          upr = predict(mod_SH_time1_SPLI_Simp2, newdata = data_pooled ,se.fit = T,type = "response" )$fit +
                            predict(mod_SH_time1_SPLI_Simp2, newdata = data_pooled,se.fit = T,type = "response" )$se.fit ,
                          lwr = predict(mod_SH_time1_SPLI_Simp2, newdata = data_pooled ,se.fit = T ,type = "response")$fit -
                            predict(mod_SH_time1_SPLI_Simp2, newdata = data_pooled ,se.fit = T ,type = "response")$se.fit , 
                          data_pooled)


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

###fit a model just with the significant interaction of TEMP_REGIME with time to 
##plot the differences
modTemp<-glm(Sh_div~ splines::ns(NumDays,3)+Temp_Regime + 
               + splines::ns(NumDays,3):Temp_Regime,
             data =data_pooled, family = "gaussian")

gg.df_TEMP <- data.frame(preds = predict(modTemp, newdata = data_pooled,type = "response" ),
                         upr = predict(modTemp, newdata = data_pooled ,se.fit = T,type = "response" )$fit +
                           predict(modTemp, newdata = data_pooled,se.fit = T,type = "response" )$se.fit ,
                         lwr = predict(modTemp, newdata = data_pooled ,se.fit = T ,type = "response")$fit -
                           predict(modTemp, newdata = data_pooled ,se.fit = T ,type = "response")$se.fit , 
                         data_pooled)


ggplot(gg.df_TEMP, aes(x =NumDays, y= Sh_div,group = Temp_Regime,
                       col =Temp_Regime, fill = Temp_Regime))+
  geom_point(alpha = 0.2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  theme_bw()



###fit a model just with the significant interaction of CORRIDORS with time to 
##plot the differences
mod_corr<-glm(Sh_div~ splines::ns(NumDays,3)+Corridor + 
               + splines::ns(NumDays,3):Corridor,
             data =data_pooled, family = "gaussian")

gg.df_Cor <- data.frame(preds = predict(mod_corr, newdata = data_pooled,type = "response" ),
                         upr = predict(mod_corr, newdata = data_pooled ,se.fit = T,type = "response" )$fit +
                           predict(mod_corr, newdata = data_pooled,se.fit = T,type = "response" )$se.fit ,
                         lwr = predict(mod_corr, newdata = data_pooled ,se.fit = T ,type = "response")$fit -
                           predict(mod_corr, newdata = data_pooled ,se.fit = T ,type = "response")$se.fit , 
                         data_pooled)


ggplot(gg.df_Cor, aes(x =NumDays, y= Sh_div,group = Corridor,
                       col =Corridor, fill = Corridor))+
  geom_point(alpha = 0.2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  theme_bw()

#----------------------------------------------------------------------------------------
## ANALYSIS of Shannon diversity VARIANCE THROUGH TIME ##
#----------------------------------------------------------------------------------------

data_pooled_Sh_div_VAR<-data_pooled %>%
  dplyr::group_by(Date,Temp_Regime, Corridor, NumDays) %>% # for each group, estimate:
  dplyr::summarise(median.Sh_div = median(Sh_div), # median Shannon diversity
                   sd.val = sd(Sh_div), # Shannon diversity standard deviation
                   Sh_variance = var(Sh_div), # Shannon diversity variance
                   upr = median.Sh_div + sd.val, # Shannon diversity upper confidence interval
                   lwr = median.Sh_div - sd.val) %>%  # Shannon diversity lower confidence interval
  filter(NumDays != 0) #avoid considering the first day in which we have variance zero

####Build the the model with cubic Spline for the time variable 
##model with all interactions with the time factor (in this case NumDays=number of days of the experiment)
mod_SH_VAR_SPL1<-glm(Sh_variance~ splines::ns(NumDays,3) +
                       Temp_Regime + 
                       Corridor + 
                       splines::ns(NumDays,3):Temp_Regime + splines::ns(NumDays,3):Corridor + Temp_Regime:Corridor+
                       splines::ns(NumDays,3):Temp_Regime:Corridor, data =data_pooled_Sh_div_VAR,
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
               data =data_pooled_Sh_div_VAR)
summary(anova_var)

###extract pairwise comparison and create a dataframe with just the significant p values of comparisons
sig_comp<-as.data.frame(TukeyHSD(anova_var)$"Temp_Regime:Corridor")
sig_comp<-sig_comp %>% filter(`p adj` < 0.05)
sig_comp

##predict values from the model and plot it, to see the fit of the model to real data
gg.df_VAR <- data.frame(preds = predict(mod_SH_VAR_SPL1, newdata = data_pooled_Sh_div_VAR,type = "response" ),
                        upr = predict(mod_SH_VAR_SPL1, newdata = data_pooled_Sh_div_VAR ,se.fit = T,type = "response" )$fit +
                          predict(mod_SH_VAR_SPL1, newdata = data_pooled_Sh_div_VAR,se.fit = T,type = "response" )$se.fit ,
                        lwr = predict(mod_SH_VAR_SPL1, newdata = data_pooled_Sh_div_VAR ,se.fit = T ,type = "response")$fit -
                          predict(mod_SH_VAR_SPL1, newdata = data_pooled_Sh_div_VAR ,se.fit = T ,type = "response")$se.fit , 
                        data_pooled_Sh_div_VAR)

ggplot(gg.df_VAR, aes(x =NumDays, y= Sh_variance,group = Temp_Regime, col =Temp_Regime, fill = Temp_Regime))+
  geom_point(alpha = 0.2)+
  geom_line(aes(y=preds))+
  geom_ribbon(aes(ymin = lwr,
                  ymax = upr),colour="transparent", alpha=0.2)+
  facet_grid(~Corridor)+
  theme_bw()

ggplot(data_pooled_Sh_div_VAR, aes(x = Temp_Regime, y = Sh_variance, col = Corridor))+
  geom_boxplot()+
  theme_bw()


##code the simpler model with splines
mod_SH_VAR_SPL2<-glm(Sh_variance~ splines::ns(NumDays,3) +
                       Temp_Regime + 
                       Corridor + 
                       splines::ns(NumDays,3):Temp_Regime + splines::ns(NumDays,3):Corridor + Temp_Regime:Corridor,
                     data =data_pooled_Sh_div_VAR,
                     family = gaussian(link = "identity"))


##Likelihood Ratio Test to compare if our model selection is significant 
lmtest::lrtest(mod_SH_VAR_SPL2,mod_SH_VAR_SPL1)

