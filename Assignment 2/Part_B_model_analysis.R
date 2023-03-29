rm(list=ls())
library(GGally)
library(corrplot)
library(car)
library(reshape)
library(tidyverse)
library(magrittr)
library(dplyr)
library(betareg)
library(statmod)
library(jtools)
library(broom.mixed)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/")
}

earinfect <- read.table("Assignment 2/earinfect.txt", header = T)



#### Earinfections ####
#1)
library(car)
head(earinfect)
earinfect$swimmer <- factor(earinfect$swimmer)
earinfect$location <- factor(earinfect$location)
earinfect$age <- factor(earinfect$age, ordered = T)
earinfect$sex <- factor(earinfect$sex)


# Poisson model -----------------------------------------------------------

fit.pois <- glm(infections ~ offset(persons) + age * sex * location * swimmer, data = earinfect, family = poisson(link = 'log'))
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-age:sex:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-age:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-age:sex:location:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-age:location:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-sex:location:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-sex:location)
Anova(fit.pois, type = "III", test = "LR")
anova(fit.pois, test = "Chisq")


# Log-normal model (Gaussian) ----------------------------------------------------------

X0 <- model.matrix(~age + sex + location + swimmer + persons, data = earinfect)
head(X0)

#looking at log-normal model: 
fit.lognorm <- lm(infections ~ age + sex + location + swimmer + persons + 
                I(location:swimmer)+I(sex:location)+I(sex:age)+I(sex:swimmer)+I(age:swimmer)+I(location:swimmer:sex), 
                data = earinfect)
anova(fit.lognorm)
fit.lognorm <- update(fit.lognorm, .~.-I(age:swimmer))
anova(fit.lognorm)
fit.lognorm <- update(fit.lognorm, .~.-I(sex:swimmer))
anova(fit.lognorm)
fit.lognorm <- update(fit.lognorm, .~.-I(location:swimmer:sex))
anova(fit.lognorm)
fit.lognorm <- update(fit.lognorm, .~.-I(sex:age))
anova(fit.lognorm)
fit.lognorm <- update(fit.lognorm, .~.-I(location:swimmer))
anova(fit.lognorm)


# Plots -------------------------------------------------------------------

1-pchisq(fit.pois$deviance, df = fit.pois$df.residual, lower.tail = F) #model is appropriate
par(mfrow = c(2,2))
plot(fit.pois)

1-pchisq(fit.lognorm$deviance, df = fit.lognorm$df.residual, lower.tail = F) #model is appropriate
plot(fit.lognorm)




# Summaries ---------------------------------------------------------------

summary(fit.pois)
summary(fit.lognorm)
summary(fit.binom)

par(mfrow=c(1,2))
acf(residuals(fit.pois, type = "pearson"))
pacf(residuals(fit.pois, type = "pearson"))
# 
# # Log-distribution of persons
# earinfect$logPersons <- log(earinfect$persons)
# earinfect$loginfections <- log(earinfect$infections)
# ggpairs(earinfect)
