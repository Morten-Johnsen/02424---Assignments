rm(list=ls())
if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 2/")
  figpath <- "/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 2/figs/"
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 2")
  figpath <- "~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 2/figs/"
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 2/")
  figpath <- "C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 2/"
}
library(data.table)
library(ggplot2)
library(GGally)
library(car)
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


predf1 <- predict(fit1,type="response")
plot(data$clo)
points(predf1, col='red')
predf0 <- predict(fit0,type="response")
plot(data$clo)
points(predf0, col='red')


# Use betaregression
library(betareg)
betafit0 <- betareg(clo ~ tOut + tInOp + sex, data = data, link = "logit")
summary(betafit0)
plot(betafit0)
res <- betafit0$residuals
qqnorm(res)
qqline(res)

# betafit0 seems good

plot(data$tOut,res)
plot(data$tIn,res)
plot(factor(data$sex),res)

betafitALL <- betareg(clo ~ tOut + tInOp + sex + I(tOut^2) + I(tInOp^2) + I(tOut*tInOp) + I(tOut^(-2)), data = data, link = "logit")
summary(betafitALL)
Anova(betafitALL,type = "III")

betafit1 <- betareg(clo ~ tOut + tInOp + sex + I(tOut^2), data = data, link = "logit")
summary(betafit1)


betafit2 <- betareg(clo ~ tOut + tInOp + sex + I(tInOp^2) + I(tOut*tInOp), data = data, link = "logit")
summary(betafit2)
Anova(betafit2,type = "III")

betafit3 <- betareg(clo ~ tOut + tInOp + sex + I(tInOp^2) , data = data, link = "logit")
summary(betafit3)
Anova(betafit3,type = "III")

betafit4 <- betareg(clo ~ tOut + sex + I(tInOp^2) , data = data, link = "logit")
summary(betafit4)
Anova(betafit4,type = "III")
par(mfrow=c(2,3))
plot(betafit4)
res <- betafit4$residuals
qqnorm(res)
qqline(res)

# betafit4 seems good

AIC(betafit4)
AIC(betafit0)


# Plot predictions
pred4 <- predict(betafit4,type="response")

data$clo
par(mfrow=c(1,1))
plot(data$clo)
points(pred4, col='red')

# Not overwhelming:((((

max(pred4)

pred0 <- predict(betafit0,type="response")

data$clo
par(mfrow=c(1,1))
plot(data$clo)
points(pred0, col='red')


