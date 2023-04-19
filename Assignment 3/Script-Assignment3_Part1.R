rm(list=ls())

library(reshape)
library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(nlme)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 3/")
  figpath <- "/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 3/figs/"
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 3")
  figpath <- "~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 3/figs/"
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/02424---Assignments/Assignment 3/")
  figpath <- "C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/02424---Assignments/Assignment 3/figs"
}

clothing <- read.csv(file = "clothingFullAss03.csv", header = T)
head(clothing)

unique(clothing$day)
unique(clothing$time)
c.data <- dplyr::select(clothing, -sex,-time2, -X, -subDay)
head(c.data)

# Not very pretty plot
melt(c.data,id = c('subjId','day') )%>%
  ggplot()+
  geom_point(aes(x = day, y=value,colour = subjId))+
  facet_wrap(~variable, scales = "free")+
  theme_bw()

mean.data <- c.data %>%
  group_by(subjId) %>%
  summarise(meanClo = mean(clo))
# Okay plot :)))) 
c.data
melt(c.data,id = c('subjId','clo'))%>%
  ggplot()+
  geom_point(aes(x = value, y=clo,colour = subjId), show.legend = FALSE)+
  facet_wrap(~variable, scales = "free")+
  theme_bw()

#Another plot looking at the 10 first subjects, clearly showing a difference between subjects depending on temperature
c.data
melt(c.data,id = c('subjId','clo','day'))%>%
  filter(day == 1 & subjId %in% unique(c.data$subjId)[1:10]) %>%
  ggplot()+
  geom_point(aes(x = value, y=clo, colour = factor(subjId)))+
  geom_line(aes(x = value, y=clo, colour = factor(subjId), linetype = factor(day)), show.legend = FALSE)+
  facet_wrap(~variable, scales = "free")+
  theme_bw()+
  labs(colour = "SubjectID:", y = "Clothing insulation level", x = "Value")+
  theme(legend.position = "top")+
  ggtitle("Between subject differences in Clothing insulation level")

# Make factors
c.data$subjId <- factor(c.data$subjId)
c.data$day <- factor(c.data$day)
str(c.data)

#### ------------- Part 1: Fitting a linear mixed effect model ---------------

## Using subjId as random effect ##
# I use the following method: first fit full model and reduce. Then fit another
# model, starting with a simple model and extending. Compare the two models
# Using chisq-testing: it does not take the random effects into account
# therefore, we compare the log-likelihoods to see if there is an improvement.
# Jan approved of comparing models during reduction/extension using chisq, and 
# then testing for the model improvement with anova. In the end, we should also 
# consider whether the random effects are significant or not.


## Simple lme: use ML method so we can compare models!
fit.mm<-lme(clo~tOut+tInOp+time+day+tOut*tInOp*time*day, 
            random = ~1|subjId, data=c.data, method="ML")

fit.mm
anova(fit.mm)
drop1(fit.mm, test = "Chisq") 

# Drop tOut:tInOp:time:day
fit.mm1 <- update(fit.mm, .~.-tOut:tInOp:time:day)
anova(fit.mm1)
drop1(fit.mm1, test = "Chisq") 
logLik(fit.mm)
logLik(fit.mm1)
anova(fit.mm,fit.mm1)
# log likelihood is lower, but we continue. p-value is big :))

# drop tInOp:time:day
fit.mm2 <- update(fit.mm1, .~.-tInOp:time:day)
anova(fit.mm2)
drop1(fit.mm2, test = "Chisq")
logLik(fit.mm)
logLik(fit.mm2)
anova(fit.mm1,fit.mm2)
# log likelihood is lower, but we continue. p-value is big :))

# drop tInOp:time:day
fit.mm3 <- update(fit.mm2, .~.-tOut:tInOp:time)
anova(fit.mm3)
drop1(fit.mm3, test = "Chisq")
logLik(fit.mm)
logLik(fit.mm3)
anova(fit.mm2,fit.mm3)
# log likelihood is lower, but we continue. p-value is big :))

# drop tInOp:time:day
fit.mm4 <- update(fit.mm3, .~.-tInOp:time)
anova(fit.mm4)
drop1(fit.mm4, test = "Chisq")
logLik(fit.mm)
logLik(fit.mm4)
anova(fit.mm3,fit.mm4)
# log likelihood is lower, but we continue. p-value is big :))

# Final model!
anova(fit.mm4)

# Now use REML so we can report model parameters
fit.mm4$terms
fit.mmREML<-lme(clo ~ tOut + tInOp + time + day + tOut:tInOp + tOut:time + tOut:day + 
                  tInOp:day + time:day + tOut:tInOp:day + tOut:time:day, 
                random = ~1|subjId, data=c.data, method="REML")
anova(fit.mmREML)

## Forward instead
fit.mmfw<-lme(clo~tOut,
            random = ~1|subjId, data=c.data, method="ML")
add1(object = fit.mmfw, scope = ~.+tInOp+time+day, test = "Chisq")
# add day
fit.mmfw1<-update(fit.mmfw, .~.+day)
anova(fit.mmfw1,fit.mmfw) # new model is better

add1(object = fit.mmfw1, scope = ~.+tInOp+time+day*tOut, test = "Chisq")
# add tOut:day
fit.mmfw2<-update(fit.mmfw1, .~.+tOut:day)
anova(fit.mmfw2,fit.mmfw1) # new model is better

add1(object = fit.mmfw2, scope = ~.+tInOp+time, test = "Chisq")
# add time
fit.mmfw3<-update(fit.mmfw2, .~.+time)
anova(fit.mmfw3,fit.mmfw2) # new model is better

add1(object = fit.mmfw3, scope = ~.+tInOp+time*tOut*day, test = "Chisq")
# add day:time
fit.mmfw4<-update(fit.mmfw3, .~.+day:time)
anova(fit.mmfw4,fit.mmfw3) # new model is better

add1(object = fit.mmfw4, scope = ~.+tInOp+time*tOut*day, test = "Chisq")
# stop here

# final model:
fit.mmfw4$terms
# compare with backward:
anova(fit.mmfw4,fit.mm4)
# The forward model is better and simpler

fit.mmfwREML<-lme(clo ~ tOut + day + time + tOut:day + day:time, 
                random = ~1|subjId, data=c.data, method="REML")
anova(fit.mmfwREML) # This test does also not take random effects into account.

