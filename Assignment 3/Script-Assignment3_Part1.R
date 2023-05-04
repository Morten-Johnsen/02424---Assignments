rm(list=ls())

library(reshape)
library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(nlme)
library(gridExtra)
library(grid)
library(lattice)

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


# Make factors
c.data$subjId <- factor(c.data$subjId)
c.data$day <- factor(c.data$day)
c.data$subDay <- factor(c.data$subDay)
c.data$sex <- factor(c.data$sex)
str(c.data)


# Summarising tables
table(c.data$sex)
table(c.data$day)

#### ------------- Part 1: Fitting a linear mixed effect model ---------------

## Using subjId as random effect ##
# I use the following method: first fit full model and reduce. Then fit another
# model, starting with a simple model and extending. Compare the two models
# Using chisq-testing: it does not take the random effects into account
# therefore, we compare the log-likelihoods to see if there is an improvement.
# Jan approved of comparing models during reduction/extension using chisq, and 
# then testing for the model improvement with anova. In the end, we should also 
# consider whether the random effects are significant or not.


# make models with individual slopes. Choose the best one.
# only use forward selection.

## Simple lme: use ML method so we can compare models!
fit.mmfw0<-lme(clo~1,
              random = ~1|subjId, data=c.data, method="ML")
add1(object = fit.mmfw0, scope = ~.+tOut+ tInOp+time+day+sex, test = "Chisq")
## Forward selection
fit.mmfw<-lme(clo~tOut,
              random = ~1|subjId, data=c.data, method="ML")
anova(fit.mmfw,fit.mmfw0)
add1(object = fit.mmfw, scope = ~.+tInOp+time+day+sex, test = "Chisq")
# add day
fit.mmfw1<-update(fit.mmfw, .~.+day)
anova(fit.mmfw1,fit.mmfw) # new model is better

add1(object = fit.mmfw1, scope = ~.+tInOp+time+sex+day*tOut, test = "Chisq")
# add tOut:day
fit.mmfw2<-update(fit.mmfw1, .~.+tOut:day)
anova(fit.mmfw2,fit.mmfw1) # new model is better

add1(object = fit.mmfw2, scope = ~.+tInOp+time+sex, test = "Chisq")
# add sex
fit.mmfw3<-update(fit.mmfw2, .~.+sex)
anova(fit.mmfw3,fit.mmfw2) # new model is better

add1(object = fit.mmfw3, scope = ~.+tInOp+time, test = "Chisq")
# add time
fit.mmfw4<-update(fit.mmfw3, .~.+time)
anova(fit.mmfw4,fit.mmfw3) # new model is better

add1(object = fit.mmfw4, scope = ~.+tInOp+day*time*tOut, test = "Chisq")
# add tInOp
fit.mmfw5<-update(fit.mmfw4, .~.+tInOp)
anova(fit.mmfw5,fit.mmfw4) # new model is better

add1(object = fit.mmfw5, scope = ~.+tInOp*day*time*tOut, test = "Chisq")

# stop here
# final model:
fit.mmfw5$terms
fit.mmfw5reml<-lme(clo ~ tOut + day + sex + time + tInOp + tOut:day,
               random = ~1|subjId, data=c.data, method="REML")
summary(fit.mmfw5reml)
tmp <- round(ranef(fit.mmfw5reml),3)
tmp$`(Intercept)`


# Make fit with individual slopes
fit.slope1 <- lme(clo ~ tOut + day + sex + time + tInOp + tOut:day, 
                    random = ~1+tOut|subjId, data=c.data, method="ML")

fit.slope2 <- lme(clo ~ tOut + time + day + sex + tOut:day, 
                    random = ~1+time|subjId, data=c.data, method="ML")

fit.slope3 <- lme(clo ~ tOut + time + day + sex + tOut:day, 
                    random = ~1+day|subjId, data=c.data, method="ML")

fit.slope4 <- lme(clo ~ tOut + time + day + sex + tOut:day, 
                  random = ~1+sex|subjId, data=c.data, method="ML")

fit.slope5 <- lme(clo ~ tOut + time + day + sex + tOut:day, 
                  random = ~1+tInOp|subjId, data=c.data, method="ML")

logLik(fit.slope1)
logLik(fit.slope2)
logLik(fit.slope3)
logLik(fit.slope4)
logLik(fit.slope5)
# fit.slope3 wins. 
# fit.slope1 is second. 

# We build a model that has a slope on the variable day too. 
# Build from scratch using forward selection
fit.slopefw1 <- lme(clo ~ tOut, 
                  random = ~1+day|subjId, data=c.data, method="ML")
# Won't converge with day as slope!
# Use tOut instead
fit.slopefw1 <- lme(clo ~ 1, 
                    random = ~1+tOut|subjId, data=c.data, method="ML")

add1(object = fit.slopefw1, scope = ~.+tOut+tInOp+time+day+sex, test = "Chisq")
# add day
fit.slopefw2<-update(fit.slopefw1, .~.+day)
anova(fit.slopefw1,fit.slopefw2) # new model is better

add1(object = fit.slopefw2, scope = ~.+tOut+tInOp+time+sex, test = "Chisq")
# add tOut
fit.slopefw3<-update(fit.slopefw2, .~.+tOut)
anova(fit.slopefw2,fit.slopefw3) # new model is better

add1(object = fit.slopefw3, scope = ~.+tInOp+time+sex+tOut*day, test = "Chisq")
# add tOut:day
fit.slopefw4<-update(fit.slopefw3, .~.+tOut:day)
anova(fit.slopefw3,fit.slopefw4) # new model is better

add1(object = fit.slopefw4, scope = ~.+tInOp+time+sex, test = "Chisq")
# add time
fit.slopefw5<-update(fit.slopefw4, .~.+time)
anova(fit.slopefw4,fit.slopefw5) # new model is better

add1(object = fit.slopefw5, scope = ~.+tInOp+sex+time*day*tOut, test = "Chisq")
# add sex
fit.slopefw6<-update(fit.slopefw5, .~.+sex)
anova(fit.slopefw5,fit.slopefw6) # new model is better

add1(object = fit.slopefw6, scope = ~.+tInOp+time*day*tOut, test = "Chisq")

# Stop here.
# Final model with slope on random effects is fit.slopefw6. 
# Compare with the other model:
anova(fit.mmfw5,fit.slopefw6)
# Model with slope is better!
fit.slopefw6$terms


final.model.ML1 <- fit.slopefw6

# Fit REML model:
final.model.REML1<-lme(clo ~ day + tOut + time + sex + day:tOut, 
                       random = ~1+tOut|subjId, data=c.data, method="REML")
summary(final.model.REML1)
# Look at random effects:
round(ranef(final.model.REML1),3)
# We can see that there is a big difference between the effect of tOut on 
# different subjects!


#   -----------------------------------------------------------------------
######### Fit a mixed effect model that include subjId and day #########
# Effects needs to be related in order to be nested
# In this case we should nest to subjId

# Use same term as fit.mm3$terms as beginning:
# + tInOp + time + day + tOut:tInOp + tOut:time + tOut:day + tInOp:day + time:day + tOut:tInOp:day + tOut:time:day

# Forward selectiong:
fit.mm.nest <- lme(clo ~ 1,
                   random = ~1|subjId/day,        
                   data = c.data, method = "ML")

add1(object = fit.mm.nest, scope = ~.+ tOut + tInOp + time + sex, test = "Chisq") # ? subDay
# add tInOp
fit.mm.nest1 <- update(fit.mm.nest, .~.+ tInOp)
anova(fit.mm.nest1, fit.mm.nest) # new model is better

add1(object = fit.mm.nest1, scope = ~.+ tOut + time + sex, test = "Chisq") # ? subDay
# add sex 
fit.mm.nest2 <- update(fit.mm.nest1, .~.+ sex)
anova(fit.mm.nest2, fit.mm.nest1) # new model is better

add1(object = fit.mm.nest2, scope = ~.+ time + tOut + tInOp*sex, test = "Chisq") # ? subDay
# add tInOp:sex
fit.mm.nest3 <- update(fit.mm.nest2, .~.+ tInOp:sex)
anova(fit.mm.nest3, fit.mm.nest2) # new model is better

add1(object = fit.mm.nest3, scope = ~.+time + tOut , test = "Chisq")
# Stop here!

# final model:
fit.mm.nest3$terms

# The forward model is better and simpler
fit.mm.nest.fwREML<-lme(clo ~ tInOp + sex + tInOp:sex, 
                        random =  ~1|subjId/day, 
                        data = c.data, method = "REML")
ranef(fit.mm.nest.fwREML)
anova(fit.mm.nest.fwREML) # This test does also not take random effects into account.
summary(fit.mm.nest.fwREML)


# Final models for later comparison
final.model.REML2 <- fit.mm.nest.fwREML
final.model.ML2 <- fit.mm.nest3

# Compare with model from 1.2:

anova(final.model.ML2,final.model.ML1)
# Conclusion: final.model.ML2 is better

# Visualise ---------------------------------------------------------------

plot(final.model.ML2)

coef(final.model.ML2)

# Precit clothing insulation levels
final.model.ML2$fit <- predict(final.model.ML2)

# Residuals versus sex
plot(final.model.ML2, resid(., type = "p") ~ fitted(.) | sex, abline = 0)

# clothing insultaion level facetted on sex
plot(final.model.ML2, clo ~ fitted(.) | sex, abline = c(0,1))


### Todo:
# Ude subDay as approximation


#   -----------------------------------------------------------------------
######### Fit a mixed effect model with within day auto-correlation #########
# subjId is random effect
# Only consider random intercept - not slope AND intercept
# Remember to validate which correlation term should be used. 
# Confirm with AIC or BIC

# Forward selection:
# insists that the grouping variable for the random effects and for the correlation be the same
fit.mm.autocor <- lme(clo ~ 1,
                      random = ~1|subDay,
                      correlation = corAR1(form = ~time|subDay),
                      data = c.data, method = "ML")

add1(object = fit.mm.autocor, scope = ~.+ tOut + tInOp + time + sex, test = "Chisq")
# add tOut
fit.mm.autocor1 <- update(fit.mm.autocor, .~.+ tOut)
anova(fit.mm.autocor1, fit.mm.autocor) # new model is better

add1(object = fit.mm.autocor1, scope = ~.+ tInOp + time + sex, test = "Chisq")
# add sex
fit.mm.autocor2 <- update(fit.mm.autocor1, .~.+ sex)
anova(fit.mm.autocor2, fit.mm.autocor1) # new model is better

add1(object = fit.mm.autocor2, scope = ~.+  time + tInOp + tOut*tInOp*time*sex, test = "Chisq")
# add tOut*sex
fit.mm.autocor3 <- update(fit.mm.autocor2, .~.+ tOut*sex) # + sex*tInOp
anova(fit.mm.autocor3, fit.mm.autocor2) # new model is better

add1(object = fit.mm.autocor3, scope = ~.+ tInOp + time + tOut*tInOp*time*sex + I(tOut^2) + I(tInOp^2) + I(time^2), test = "Chisq")
# Nothing more to add


# final model:
fit.mm.autocor3$terms

# The forward model is better and simpler
fit.mm.autocor.fwREML <- lme(clo ~ tOut + sex + tOut*sex,
                      random = ~1|subDay,
                      correlation = corAR1(form = ~time|subDay),
                      data = c.data, method = "REML")
ranef(fit.mm.autocor.fwREML)
anova(fit.mm.autocor.fwREML) # This test does also not take random effects into account.
summary(fit.mm.autocor.fwREML)

# Final models for later comparison
final.model.REML3 <- fit.mm.autocor.fwREML
final.model.ML3 <- fit.mm.autocor3

# For interpretation:
subID_interp=12
ranef.data = data.frame(ranef(fit.mm.autocor.fwREML),subDay = 0:135)
c.data.ranef = merge(c.data,ranef.data,by = "subDay")

c.data.ranef[c.data.ranef$subjId == subID_interp,]

samefit_wo_AR <- lme(clo ~ tOut + sex + tOut*sex,
                             random = ~1|subDay,
                             data = c.data, method = "REML")
ranef.data_wo_AR = data.frame(ranef(samefit_wo_AR),subDay = 0:135)
c.data.ranef_wo_AR = merge(c.data.ranef,ranef.data_wo_AR,by = "subDay")

c.data.ranef_wo_AR[c.data.ranef_wo_AR$subjId == subID_interp,]

samefit_wo_ARML <- lme(clo ~ tOut + sex + tOut*sex,
                     random = ~1|subDay,
                     data = c.data, method = "ML")
anova(final.model.ML3,samefit_wo_ARML)


# COMPARE ALL THREE MODELS: -----------------------------------------------------
anova(final.model.ML1,final.model.ML2)
anova(final.model.ML2,final.model.ML3)

# forestmodel::forest_model(final.model.ML1, final.model.ML2, final.model.ML3)
AIC(final.model.ML1)
AIC(final.model.ML2)
AIC(final.model.ML3)
# AIC is lowest for the last model with autocorrelation


# Visualise ---------------------------------------------------------------

plot(final.model.ML3)

coef(final.model.ML3)

# Precit clothing insulation levels
final.model.ML3$fit <- predict(final.model.ML3)

# Residuals versus sex
plot(final.model.ML3, resid(., type = "p") ~ fitted(.) | sex, abline = 0)

# clothing insultaion level facetted on sex
plot(final.model.ML3, clo ~ fitted(.) | sex, abline = c(0,1))

# This also shows higher variance of fitted values vs. actual values for females (than male)

par(mfrow = c(1, 1))
acf(residuals(final.model.ML3, retype="normalized"), main = "Auto-correlation for final.model.ML3 (ACF)", lag.max = 6)
pacf(residuals(final.model.ML3, retype="normalized"), main = "Auto-correlation for final.model.ML3 (PACF)", lag.max = 6)


# QQplot
c.data$res <- resid(final.model.ML3, type = "p")

qqplot <- ggplot(c.data, aes(sample = res, colour = factor(sex))) +
  stat_qq() +
  stat_qq_line() + 
  theme_minimal()+
  labs(title = "Normal QQ", subtitle = "Within-day autocorrelation",
       x = "Theoretical quantiles", y = "Std. Persons residuals", color = "Sex") 
ggsave(filename = file.path(figpath, "qqplot.png"), plot = qqplot, height = 5, width = 7.5)



# Predicted versus actual -------------------------------------------------
ML3 <- ggplot(c.data, aes(x = predict(final.model.ML3), y = clo)) +
  geom_smooth(method = "lm") +
  geom_point() +
  labs(x='Predicted Values', y='Actual Values', title='Predicted vs. Actual Values', subtitle = "With within day auto-correlation") +
  theme_minimal() 
ML2 <- ggplot(c.data, aes(x = predict(final.model.ML2), y = clo)) +
  geom_smooth(method = "lm") +
  geom_point() +
  labs(x='Predicted Values', y='Actual Values', title='Predicted vs. Actual Values', subtitle = "With subjId and day as random effects") +
  theme_minimal() 
ML1 <- ggplot(c.data, aes(x = predict(final.model.ML1), y = clo)) +
  geom_smooth(method = "lm") +
  geom_point() +
  labs(x='Predicted Values', y='Actual Values', title='Predicted vs. Actual Values', subtitle = "With subjId as random effects") +
  theme_minimal() 

ggsave(filename = file.path(figpath, "predictVSactual1.png"), plot = ML1, height = 5, width = 7.5)
ggsave(filename = file.path(figpath, "predictVSactual2.png"), plot = ML2, height = 5, width = 7.5)
ggsave(filename = file.path(figpath, "predictVSactual3.png"), plot = ML3, height = 5, width = 7.5)



# With pred. intervals ----------------------------------------------------

newdat <- data.frame( tOut = seq(min(c.data$tOut), max(c.data$tOut), length = nrow(c.data)), 
                      tInOp = seq(min(c.data$tInOp), max(c.data$tInOp), length = nrow(c.data)),
                      sex = c(rep("male", 400), rep("female", 403)),
                      subjId = sort(c(rep(seq(48, 100, 10), 133), 48, 58, 68, 78, 98)),
                      time = seq(min(c.data$time), max(c.data$time), length = nrow(c.data)),
                      day = sort(c(rep(c(1, 2, 3), 267), 1, 1)),
                      subDay = sort(c(rep(c(1, 2, 3), 267), 1, 1)))  

# Make factors
newdat$subjId <- factor(newdat$subjId)
newdat$day <- factor(newdat$day)
newdat$subDay <- factor(newdat$subDay)
newdat$sex <- factor(newdat$sex)

newdat$pred <- predict(final.model.ML3, newdat, level = 0)

Designmat <- model.matrix(eval(eval(final.model.ML3$call$fixed)[-2]), newdat[-ncol(newdat)])

# Compute standard error for predictions
predvar <- diag(Designmat %*% final.model.ML3$varFix %*% t(Designmat))
newdat$SE <- sqrt(predvar) 
newdat$SE2 <- sqrt(predvar+final.model.ML3$sigma^2) # sigma = stdDev residual = SE*2 gives 95% confidence intervals on predictions
newdat$lowerCI<-newdat$pred-(2*newdat$SE2)
newdat$upperCI<-newdat$pred+(2*newdat$SE2)


c.data$pred <- predict(final.model.ML3)
ML3.pred <- ggplot(c.data,aes(x=clo,y=pred)) + 
  geom_point() +
  geom_ribbon(aes(ymin=pred-2*newdat$SE,ymax=pred+2*newdat$SE),alpha=0.3,fill="blue") +
  geom_ribbon(aes(ymin=pred-2*newdat$SE2,ymax=pred+2*newdat$SE2),alpha=0.2,fill="red") +
  labs(x='Actual Values', y='Predicted Values', title='Predicted vs. Actual Values of Clothing Insulation Level',
       subtitle = "Within day auto-correlation") +
  theme_minimal()
ggsave(filename = file.path(figpath, "predictVSactual3_1.png"), plot = ML3.pred, height = 5, width = 7.5)


# Facet on sex
ML3.pred.sex <- ggplot(c.data,aes(x=clo,y=pred)) + 
  geom_point() +
  geom_ribbon(aes(ymin=pred-2*newdat$SE,ymax=pred+2*newdat$SE),alpha=0.3,fill="blue") +
  geom_ribbon(aes(ymin=pred-2*newdat$SE2,ymax=pred+2*newdat$SE2),alpha=0.2,fill="red") +
  labs(x='Actual Values', y='Predicted Values', title='Predicted vs. Actual Values of Clothing Insulation Level',
       subtitle = "Within day auto-correlation") +
  theme_minimal() +
  facet_wrap(vars(sex))
ggsave(filename = file.path(figpath, "predictVSactual3_1_sex.png"), plot = ML3.pred.sex, height = 5, width = 7.5)

boxplot(c.data$pred[c.data$sex == "male"], c.data$clo[c.data$sex == "male"])




### How much does the random effects explain the left over variance?
## ML1 (tOut):
0.01271465/(0.01271465 + 0.08153435) * 100
# > 13.5 %

