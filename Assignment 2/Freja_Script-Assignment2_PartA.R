rm(list=ls())
if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 2/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 2")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 2/")
}
library(data.table)
library(ggplot2)
library(GGally)
#load data
data <- data.table(read.csv("clothing.csv",header=TRUE))

# Look at data
head(data)
tail(data)
summary(data)

# Remove index column
data$X <- NULL


ggpairs(data)
par(mfrow=c(2,3))
plot(clo~tOut,data=data,main="tOut")
plot(clo~tInOp,data=data,main="tInOp")
boxplot(clo~sex,data=data,main="sex")
boxplot(clo~subjId,data=data,main="subjId")
plot(clo~time,data=data,main="time")
boxplot(clo~day,data=data,main="day")

# Based on the plot:
#   - There seems to be a correlation between temperature and clothing
#   - Also seems to be effect of sex
#   - No effect from day or time or subjId

# Question: Does it make sense to use subjid???

# Let's look at y before choosing distribution
data$clo
# Values between 0 and 1 - kind of like percentage?
#     - Can't use binomial since it is counting number of successes
#     - According to the internet, use beta distribution
# Correction: according to ?binomial, the y-vector can be given as a numerical vector
# with values between 0 and 1 (proportion of successful cases). Though, then I have to give
# the total number of cases in weights (????)

# It is wrong to set the number of trials to 100, because it is used to calculate p values.
# I try it anyways:
data$weight <- 100
fit0 <- glm(clo ~ tOut + tInOp + sex + time + day, data = data, 
            family = binomial(link = "logit"), weights = data$weight)
summary(fit0)
# All parameters are significant, which I don't believe. 

# When not specifying trials:
fit1 <- glm(clo ~ tOut + tInOp + sex + time + day, data = data, 
            family = binomial(link = "logit"))
summary(fit1)
# Not all parameters are significant, which is more reasonable

# Remove 'day'
fit2 <- glm(clo ~ tOut + tInOp + sex + time, data = data, 
            family = binomial(link = "logit"))
summary(fit2)

# Remove 'tInOp'
fit3 <- glm(clo ~ tOut + sex + time, data = data, 
            family = binomial(link = "logit"))
summary(fit3)

# Remove 'time'
fit4 <- glm(clo ~ tOut + sex, data = data, 
            family = binomial(link = "logit"))
summary(fit4)

(pval <- 1 - pchisq(67.333,800))
# Good

par(mfrow=c(2,2))
plot(fit4)
# Residuals look good

# Use betaregression
library(betareg)
betafit0 <- betareg(clo ~ tOut + tInOp + sex + time + day, data = data, link = "logit")
summary(betafit0)

# Remove day
betafit1 <- betareg(clo ~ tOut + tInOp + sex + time, data = data, link = "logit")
summary(betafit1)

plot(betafit1)
