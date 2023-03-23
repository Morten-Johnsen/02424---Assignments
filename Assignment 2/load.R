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
c.data2 <- dplyr::select(clothing, -X)
c.data2 %<>% 
  mutate(subjId = factor(subjId))
fit.gamma2 <- glm(clo ~ tOut + tInOp + subjId, data = c.data2, family = Gamma(link = "cloglog"))
anova(fit.gamma2, test = "Chisq")
#subject id is highly significant -> normally this would be an indicator to use a mixed model instead.
#the residual deviance for this model is a lot lower than for the sex-based model above.
c.data2$pearson <- residuals(fit.gamma2, type = "pearson")

ggplot(c.data2)+
  geom_boxplot(aes(x = subjId, y = pearson))+
  theme_bw()
#still not constant variance.

par(mfrow=c(2,2))
plot(fit.gamma2)

#### 5) within day autocorrelation
acf(residuals(fit.gamma2))
pacf(residuals(fit.gamma2))
par(mfrow=c(1,1))
#still seeing autocorrelation
c.data2 %>%
  group_by(subjId, day, time) %>%
  arrange(subjId, day, time) %>%
  summarise(autocor = acf(pearson, plot = F)$acf) -> test

hist(test$autocor)

ggplot(c.data2)+
  geom_point(aes(x = subjId, y = pearson, fill = day))

acf <- c()
lag <- c()
subject <- c()

for (i in unique(c.data2$subjId)){
  c.data2 %>%
    filter(subjId == i) -> tmp
  for (j in unique(tmp$day)){
    tmp %>%
      filter(day == j) %>%
      arrange(time) %>%
      select(pearson) %>%
      acf(plot = F) -> acf.tmp
    acf <- c(acf, acf.tmp$acf)
    lag <- c(lag, 0:(length(acf.tmp$acf)-1))
    subject <- c(subject, rep(i, length(acf.tmp$acf)))
  }
}

ggplot(data.frame("acf" = acf, "lag" = lag, "subject" = subject))+
  geom_col(aes(x = lag, y = acf, fill = subject), position = "dodge")+
  theme_bw()+
  theme(legend.position = "none")+
  scale_fill_manual(values = rep("black",length(unique(subject))))+
  geom_hline(aes(yintercept = qnorm(1-0.05/2)/sqrt(dim(c.data2)[1]), colour = "95% significance level"), linetype = "dashed", colour = "royalblue1", size = 0.8)+
  geom_hline(aes(yintercept = -qnorm(1-0.05/2)/sqrt(dim(c.data2)[1])), linetype = "dashed", colour = "royalblue1", size = 0.8)+
  scale_x_continuous(n.breaks = 6)

#6) Optimal weight/dispersion parameter
glm.gamma.w <- function(theta){
  y <- c.data2$clo
  w.male <- theta[5]
  w.female <- theta[6]
  w <- numeric(dim(c.data2)[1])
  w[c.data2$sex == "male"] <- w.male
  w[c.data2$sex == "female"] <- w.female
  
  eta <- theta[1] + theta[2] * c.data2$tOut + theta[3] * c.data2$tInOp + theta[4] * as.numeric(c.data2$sex == "male")
  mu <- 1-exp(-exp(eta))
  
  d <- 2*(y/mu - log(y/mu) - 1)
  return(1/2 * sum(w*d))
}

manual.fit <- nlminb(start = c(0,0,0,0,1,1), objective = glm.gamma.w)
manual.fit$par

#7) Profile likelihood
glm.gamma.w.pf <- function(w1,w2){
  y <- c.data2$clo
  w.male <- w1
  w.female <- w2
  w <- numeric(dim(c.data2)[1])
  w[c.data2$sex == "male"] <- w.male
  w[c.data2$sex == "female"] <- w.female
  tmp.func <- function(theta){
    eta <- theta[1] + theta[2] * c.data2$tOut + theta[3] * c.data2$tInOp + theta[4] * as.numeric(c.data2$sex == "male")
    mu <- 1-exp(-exp(eta))
    
    w <- w*mu^2
    d <- 2*(y/mu - log(y/mu) - 1)
    return(1/2 * sum(w*d))
  }
  
  fit.tmp <- nlminb(start = c(0,0,0,0), objective = tmp.func)
  return(fit.tmp$objective)
}

w1 <- seq(-5, -1, length.out = 100)
w2 <- seq(-12, -8, length.out = 100)
z <- outer(w1, w2, FUN = Vectorize(function(w1,w2) glm.gamma.w.pf(w1,w2)))

contour(w1, w2, z)

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
