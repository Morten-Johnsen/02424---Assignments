rm(list=ls())
library(GGally)
library(corrplot)
library(car)
library(reshape)
library(tidyverse)
library(dplyr)
library(betareg)
library(statmod)
library(jtools)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/")
}

earinfect <- read.table("Assignment 2/earinfect.txt", header = T)
clothing <- read.csv(file = "Assignment 2/clothing.csv", header = T)

#### Clothing ####
#1)
c.data <- dplyr::select(clothing, -subjId, -day, -time, -X)
head(c.data)

melt(c.data) %>%
  ggplot()+
  geom_histogram(aes(x = value, fill = sex), position = "identity", alpha = 0.4)+
  scale_fill_manual(values = c("blue", "orange"))+
  facet_wrap(~variable, scales = "free")+
  theme_bw()

#### 1) Distributions
beta.dist <- function(theta){
  return(-sum(dbeta(c.data$clo, shape1 = theta[1], shape2 = theta[2], log = T)))
}

norm.dist <- function(theta){
  return(-sum(dnorm(c.data$clo, mean = theta[1], sd = theta[2], log = T)))
}

gamma.dist <- function(theta){
  return(-sum(dgamma(c.data$clo, shape = theta[1], rate = theta[2], log = T)))
}

invgaus.dist <- function(theta){
  return(-sum(dinvgauss(c.data$clo, mean = theta[1], dispersion = theta[2], log =T)))
}

lnorm.dist <- function(theta){
  return(-sum(dlnorm(x=c.data$clo, meanlog = theta[1], sdlog = theta[2], log = T)))
}

beta.hat <- nlminb(start = c(1,1), objective = beta.dist)
norm.hat <- nlminb(start = c(1,1), objective = norm.dist)
gamma.hat <- nlminb(start = c(1,1), objective = gamma.dist)
invgaus.hat <- nlminb(start = c(1,1), objective = invgaus.dist)
lnorm.hat <- nlminb(start = c(1,1), objective = lnorm.dist)

ggplot(c.data)+
  geom_histogram(aes(x = clo, y = after_stat(density)), colour = "white", position = "identity", alpha = 0.4)+
  scale_fill_manual(values = c("blue", "orange"))+
  stat_function(aes(colour = "1: beta"), fun = dbeta, args = list(shape1 = beta.hat$par[1],shape2 = beta.hat$par[2]))+
  stat_function(aes(colour = "2: normal"), fun = dnorm, args = list(mean = norm.hat$par[1], sd = norm.hat$par[2]))+
  stat_function(aes(colour = "3: gamma"), fun = dgamma, args = list(shape = gamma.hat$par[1],rate = gamma.hat$par[2]))+
  stat_function(aes(colour = "4: gamma"), fun = dinvgauss, args = list(mean = invgaus.hat$par[1],dispersion = invgaus.hat$par[2]))+
  stat_function(aes(colour = "5: log-normal"), fun = dlnorm, args = list(meanlog = lnorm.hat$par[1],sdlog = lnorm.hat$par[2]))+
  theme_bw()+
  labs(y = "", colour = "Distribution")+
  scale_colour_manual(values = c("blue", "orange", "black", "red", "purple", "grey"), labels = c(paste0("Beta [AIC: ",round(2*beta.hat$objective + 2*length(beta.hat$par),3),"]"),
                                                                        paste0("Normal [AIC: ",round(2*norm.hat$objective + 2*length(norm.hat$par),3),"]"), 
                                                                        paste0("Gamma [AIC: ",round(2*gamma.hat$objective + 2*length(gamma.hat$par),3),"]"),
                                                                        paste0("Inv Gauss [AIC: ",round(2*invgaus.hat$objective + 2*length(invgaus.hat$par),3),"]"),
                                                                        paste0("Log-normal [AIC: ",round(2*lnorm.hat$objective + 2*length(lnorm.hat$par),3),"]")))+
  ggtitle("Clothing insulation level")

#gamma is best (HUSK clo må gerne være > 1, det er ikke en fraktion).
fit.gamma <- glm(clo ~ tOut + tInOp + sex, data = c.data, family = Gamma(link = "cloglog"))
#goodness of fit:
pchisq(deviance(fit.gamma), df = dim(c.data)[1] - length(coefficients(fit.gamma)), lower.tail = F)
#er passende
summary(fit.gamma)


#### 2) residual analysis
par(mfrow=c(2,2))
plot(fit.gamma)

#gender-specific residual analysis
c.data$pred <- predict(fit.gamma)
c.data$pearson <- residuals(fit.gamma, type = "pearson")

#Check whether residuals are i.i.d.
ggplot(c.data)+
  geom_boxplot(aes(x = sex, y = pearson, fill = sex))+
  theme_bw()+
  ggtitle("Residual variation difference between genders")
#=> Residuals are not identically distributed. This will be addressed later.

par(mfrow=c(1,2))
acf(c.data$pearson)
pacf(c.data$pearson)
#=> residuals are not independent.

#### 3) Model interpretation
plot_summs(fit.gamma)

#### 4) Fitting the model using subjId instead of sex
c.data2 <- dplyr::select(clothing, -sex, -day, -time, -X)
fit.gamma2 <- glm()

#### Earinfections ####
#1)
library(car)
head(earinfect)
earinfect$swimmer <- factor(earinfect$swimmer)
earinfect$location <- factor(earinfect$location)
earinfect$age <- factor(earinfect$age, ordered = T)
earinfect$sex <- factor(earinfect$sex)

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


summary(fit.pois)
#check whether the model is appropriate or if there is overdispersion.
pchisq(fit.pois$deviance, df = fit.pois$df.residual, lower.tail = F) #model is appropriate
par(mfrow = c(2,2))
plot(fit.pois)
#for poisson glm with overdispersion see p. 152

#looking at binomial model:
earinfect$resp <- cbind(earinfect$infections, earinfect$persons - earinfect$infections)
fit.logis <- glm(resp ~ age + sex + location + swimmer, data = earinfect, family = binomial(link = "logit"))
anova(fit.logis, test = "Chisq")
Anova(fit.logis, type = "III", test.statistic = "LR")
fit.logis <- update(fit.logis, .~.-sex)
anova(fit.logis, test = "Chisq")
Anova(fit.logis, type = "III", test.statistic = "LR")
fit.logis <- update(fit.logis, .~.-age)
anova(fit.logis, test = "Chisq")
Anova(fit.logis, type = "III", test.statistic = "LR")
fit.logis <- update(fit.logis, .~.-swimmer)

add1(fit.logis, scope=~.+I(location:swimmer)+I(sex:location)+I(sex:age)+I(sex:swimmer)+I(age:swimmer)+I(location:swimmer:sex), test = "Chisq")

summary(fit.logis)
pchisq(fit.logis$deviance, df = fit.logis$df.residual, lower.tail = F) #model is highly appropriate
par(mfrow=c(2,2))
plot(fit.logis)
earinfect$logis.pred <- predict(fit.logis, type = "response")
ggplot(earinfect, aes(x = as.numeric(location)))+
  geom_point(aes(y = infections/persons))+
  geom_line(aes(y = logis.pred))+
  theme_bw()+
  scale_x_continuous(breaks = c(1,2), labels = c("beach", "non beach"))+
  ggtitle("Logistic regression model")
