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
c.data <- dplyr::select(clothing, -time2, -X)
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
p1 <- melt(c.data,id = c('subjId','clo','day'))%>%
  filter(day == 1 & subjId %in% unique(c.data$subjId)[1:10]) %>%
  ggplot()+
  geom_point(aes(x = value, y=clo, colour = factor(subjId)))+
  geom_line(aes(x = value, y=clo, colour = factor(subjId), linetype = factor(day)), show.legend = FALSE)+
  facet_wrap(~variable, scales = "free")+
  theme_bw()+
  labs(colour = "SubjectID:", y = "Clothing insulation level", x = "Value")+
  theme(legend.position = "top")+
  ggtitle("Between subject differences in Clothing insulation level")
ggsave(filename = file.path(figpath, "subjIdDifferences.png"), plot = p1)


p2 <- melt(c.data,id = c('subjId','clo','day'))%>%
  filter(day == 1 & subjId %in% unique(c.data$subjId)[1:10]) %>%
  ggplot()+
  geom_histogram(aes(x = value, fill = factor(subjId))) +
  # geom_point(aes(x = value, y=clo, colour = factor(subjId)))+
  # geom_line(aes(x = value, y=clo, colour = factor(subjId), linetype = factor(day)), show.legend = FALSE)+
  facet_wrap(~subjId, scales = "free")+
  theme_bw()+
  labs(colour = "SubjectID:", y = "Clothing insulation level", x = "Value")+
  theme(legend.position = "top")+
  ggtitle("Between subject differences in Clothing insulation level")
p2


melt(c.data,id = c('subjId','clo','day'))%>%
  filter(day == 1 & subjId %in% unique(c.data$subjId)[1:10]) %>%
  ggplot()+
  geom_point(aes(x = value, y = clo, fill = factor(subjId), size = value), 
             alpha = 0.5, shape = 21, color = "black") +
  scale_size(range = c(.1, 10), name = "") +
  theme_minimal() +
  theme(legend.position = "top") +
  labs(y = "Clothing insulation level", 
       x = "Value", 
       title = "Between subject differences in Clothing insulation level") +
  guides(fill = guide_legend(title = "subjId"))



# Make factors
c.data$subjId <- factor(c.data$subjId)
c.data$day <- factor(c.data$day)
c.data$subDay <- factor(c.data$subDay)
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
logLik(fit.mm4)

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




######### Fit a mixed effect model that include subjId and day #########
# Effects needs to be related in order to be nested
# In this case we should nest to subjId


# Use same term as fit.mm4$terms as beginning:
# c.data$f <- with(c.data, subjId:day)
fit.mm.nest <- lme(clo ~ tOut + tInOp + time + day + tOut:tInOp + tOut:time + tOut:day + tInOp:day + time:day + tOut:tInOp:day + tOut:time:day,
                   random = ~1 + day|subjId,         # the effect of day can be different for each subjectID  (eller skal det v?re omvendt?)
                   # ~ 1 | f,
                   # list(subjId = ~ 1, day = ~ 1), # specifies that subjId is nested within day
                   data = c.data, method="ML") #skal det ikke være 1|subjId/day for at det er nested?


fit.mm.nest
anova(fit.mm.nest)
drop1(fit.mm.nest, test = "Chisq") 

# Drop tOut:tInOp:day
fit.mm.nest1 <- update(fit.mm.nest, .~.-tOut:tInOp:day)
anova(fit.mm.nest1)
drop1(fit.mm.nest1, test = "Chisq") 
logLik(fit.mm.nest)
logLik(fit.mm.nest1)
anova(fit.mm.nest,fit.mm.nest1)
# log likelihood is lower, but we continue. p-value is big :))

# Drop tOut:tInOp
fit.mm.nest2 <- update(fit.mm.nest1, .~.-tOut:tInOp)
anova(fit.mm.nest2)
drop1(fit.mm.nest2, test = "Chisq") 
logLik(fit.mm.nest1)
logLik(fit.mm.nest2)
anova(fit.mm.nest1,fit.mm.nest2)
# log likelihood is smaller, but we continue. p-value is big :))

# Drop tOut:time:day
fit.mm.nest3 <- update(fit.mm.nest2, .~.-tOut:time:day)
anova(fit.mm.nest3)
drop1(fit.mm.nest3, test = "Chisq") 
logLik(fit.mm.nest2)
logLik(fit.mm.nest3)
anova(fit.mm.nest2,fit.mm.nest3)
# log likelihood is lower, but we continue. p-value is big :))

# Drop tOut:tInOp:time:day
fit.mm.nest4 <- update(fit.mm.nest3, .~.-tOut:time - time:day)
anova(fit.mm.nest4)
drop1(fit.mm.nest4, test = "Chisq") 
logLik(fit.mm.nest3)
logLik(fit.mm.nest4)
anova(fit.mm.nest3,fit.mm.nest4)
# log likelihood is lower, but we continue. p-value is big :))

# Drop tOut:tInOp:time:day
fit.mm.nest5 <- update(fit.mm.nest4, .~.-tOut:day)
anova(fit.mm.nest5)
drop1(fit.mm.nest5, test = "Chisq") 
logLik(fit.mm.nest4)
logLik(fit.mm.nest5)
anova(fit.mm.nest4,fit.mm.nest5)
# log likelihood is lower, but we continue. p-value is big :))

# Drop tOut:tInOp:time:day
fit.mm.nest6 <- update(fit.mm.nest5, .~.-time)
anova(fit.mm.nest6)
drop1(fit.mm.nest6, test = "Chisq") 
logLik(fit.mm.nest5)
logLik(fit.mm.nest6)
anova(fit.mm.nest5,fit.mm.nest6)
# log likelihood is lower, but we continue. p-value is big :))


#### Part 2 ####
library(lme4)
fit0 <- lmer(clo~sex+(1|subjId),data=c.data,REML=FALSE)
summary(fit0)

X <- cbind(1, as.numeric(c.data$sex == "male"))
Z <- dummy(c.data$subjId, levelsToKeep = unique(c.data$subjId))
y <- c.data$clo
dim(Z)
#Simultaneous estimation of beta and u for known variances p. 184
beta <- solve(t(X)%*%X)%*%t(X)%*%y
beta_old <- beta
u <- matrix(0, nrow = 47)
u_old <- u

Sigma <- diag(rep(1,length(y))); Psi <- diag(rep(1,dim(Z)[2]))

iterations <- 0
#Ved ikke om det her er rigtigt - det er som om der mangler et eller andet eller at det kan gøres smartere... 
#Det virker dog til at få de rigtige parameter estimater
while ((all(abs(beta - beta_old) > 1e-9) & all(abs(u - u_old) > 1e-9)) | iterations < 1){
  
  beta_old <- beta
  u_old <- u
  iterations <- iterations + 1
  #calculate the adjusted observation
  y_adj <- y - X%*%beta
  #estimate u
  
  u <- solve(t(Z)%*%solve(Sigma)%*%Z + solve(Psi)) %*% (t(Z)%*%solve(Sigma)%*%y_adj)
  
  y_adj <- y - Z%*%u
  
  #reestimate beta
  beta <- solve(t(X)%*%solve(Sigma)%*%X)%*%t(X)%*%solve(Sigma)%*%y_adj
  
  #tmp.func <- function(theta, u, e, Z, X){
  tmp.func <- function(theta, u, e){
    Sigma <- exp(theta[1])
    Psi <- exp(theta[2])
    
    obj <- sum(dnorm(e, sd = sqrt(Sigma), log = T)) + sum(dnorm(u, sd = sqrt(Psi), log = T))
    #Sigma <- diag(rep(exp(theta[1]),length(y))); Psi <- diag(rep(exp(theta[2]),dim(Z)[2]))
    
    #V <- Sigma + Z %*% Psi %*% t(Z)
    #obj <- -0.5 * log(det(V)) - 0.5*log(det(t(X)%*%solve(V)%*%X)) - 0.5*t(e) %*% solve(V) %*% e
    return(-obj)
  }
  e <- y - X%*%beta - Z%*%u
  est <- nlminb(start = c(1,1), objective = tmp.func, u = u, e = e)
  #est <- nlminb(start = c(-2,-1), objective = tmp.func, u = u, e = e, Z = Z, X = X)
  
  Sigma <- diag(rep(exp(est$par[1]),length(y)))
  Psi <- diag(rep(exp(est$par[2]),dim(Z)[2]))

  
  if(iterations %in% c(1, 10, 100, 200, 300, 400, 600, 800, 1000, 10000)){
    cat("\nIteration: ", iterations, " done. Update difference: ",max(abs(beta-beta_old))," \n----------------------------")
  }
}
iterations
beta-beta_old
fit0
beta
cbind(ranef(fit0)$subjId, u)
sqrt(unique(diag(Psi)))
sqrt(unique(diag(Sigma)))
fit0

#Det herunder kører ikke rigtig
tmp.func <- function(theta, X, Z, y){
  
  Sigma <- exp(theta[1])
  Psi <-   exp(theta[2])
  
  beta <- matrix(theta[3:4], nrow = 2)
  u <- matrix(theta[-c(1:4)], nrow = dim(Z)[2])
  
  e <- y - X%*%beta - Z%*%u
  
  obj <- sum(dnorm(e, sd = sqrt(Sigma), log = T)) + sum(dnorm(u, sd = sqrt(Psi), log = T))
  return(-obj)
}

est_direct <- nlminb(start = c(0,0,0,0,rep(0,dim(Z)[2])), objective = tmp.func, X = X, Z = Z, y = y, 
                     control = list(trace = 1))

summary(fit0)
sqrt(exp(est_direct$par[1:2]))
est_direct$par[3:4]
cbind(ranef(fit0)$subjId,
      est_direct$par[-c(1:4)])
