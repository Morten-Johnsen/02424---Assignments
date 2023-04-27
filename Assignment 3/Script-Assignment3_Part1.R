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


#   -----------------------------------------------------------------------
# Subset data to contain oberservations for the 10 first sujects
tiny.data <- c.data %>% filter(day == 1 & subjId %in% unique(c.data$subjId)[1:10])

# Calculate mean tOut, tInop and time for each subject for day 1
tiny.mean.data <- tiny.data %>%
  group_by(subjId) %>%
  summarise(meanClo = mean(clo),
            meantOut = mean(tOut),
            meantime = mean(time),
            meantInOp = mean(tInOp))

print(tiny.mean.data)

p1 <- 
  ggplot(tiny.data, aes(x = tOut, y = clo, group = factor(subjId), color = factor(subjId)), size = 1) + 
  geom_line() +
  geom_point() +
  labs(colour = "SubjectID:", y = "Clothing insulation level", x = "Outdoor temperature")  +
  theme_minimal()
p1
p2 <- 
  ggplot(tiny.data, aes(x = tInOp, y = clo, group = factor(subjId), color = factor(subjId)), size = 1, show.legend = FALSE) + 
  geom_line() +
  geom_point() +
  labs(colour = "SubjectID:", x = "Indoor operating temperature", y = "") +
  theme_minimal() 
p2
p3 <- 
  ggplot(tiny.data, aes(x = time, y = clo, group = factor(subjId), color = factor(subjId)), size = 1, show.legend = FALSE) + 
  geom_line() +
  geom_point() +
  labs(colour = "SubjectID:", x = "Within day time measurement", y = "") +
  theme_minimal()
p3

pall <- grid.arrange(ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = TRUE, legend="bottom"),
                     top = textGrob("Between subject differences in Clothing insulation level",
                                    gp = gpar(fontsize = 12, font = 1),
                                    x = 0, hjust = 0))
ggsave(filename = file.path(figpath, "subjIdDifferences.new.png"), plot = pall, height = 5, width = 8)


#   -----------------------------------------------------------------------




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
c.data$sex <- factor(c.data$sex)
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


# TODO: include gender too
# make models with individual slopes. Choose the best one.
# only use forward selection.

## Simple lme: use ML method so we can compare models!
## Forward selection
fit.mmfw<-lme(clo~tOut,
              random = ~1|subjId, data=c.data, method="ML")
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
summary(fit.mmfw5)

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
# Look at random effects:
ranef(final.model.REML1)
# We can see that there is a big difference between the effect of tOut on 
# different subjects!


#   -----------------------------------------------------------------------
######### Fit a mixed effect model that include subjId and day #########
# Effects needs to be related in order to be nested
# In this case we should nest to subjId

# Todo: Only forward selection. Compare with final model from previously
# Tjek at det er det samme som at bruge subDay


# Use same term as fit.mm3$terms as beginning:
# + tInOp + time + day + tOut:tInOp + tOut:time + tOut:day + tInOp:day + time:day + tOut:tInOp:day + tOut:time:day

# Forward selectiong:
fit.mm.nest <- lme(clo ~ tOut,
                   random = ~1|subjId/day,        
                   data = c.data, method = "ML")

add1(object = fit.mm.nest, scope = ~.+ tInOp + time + sex, test = "Chisq") # ? subDay
# add tInOp and sex
fit.mm.nest1 <- update(fit.mm.nest, .~.+ tInOp + sex)
anova(fit.mmfw1, fit.mmfw) # new model is better

add1(object = fit.mm.nest1, scope = ~.+time + tOut + tOut*sex*time, test = "Chisq")
# add tOut:sex
fit.mm.nest2 <- update(fit.mm.nest1, .~.+tOut:sex)
anova(fit.mm.nest2, fit.mm.nest1) # new model is better

add1(object = fit.mm.nest2, scope = ~.+tInOp*tOut*sex*time, test = "Chisq")
# add tInOp:sex
fit.mm.nest3 <- update(fit.mm.nest2, .~.+tInOp:sex)
anova(fit.mm.nest3, fit.mm.nest2) # new model is better

add1(object = fit.mm.nest3, scope = ~.+tInOp*tOut*sex*time + I(tOut^2) + I(tInOp^2) + I(time^2), test = "Chisq")
# Stop here!


# final model:
fit.mm.nest3$terms

# The forward model is better and simpler
fit.mm.nest.fwREML<-lme(clo ~ tOut + tInOp + sex + tOut*sex + tInOp*sex, 
                        random =  ~1|subjId/day, 
                        data = c.data, method = "REML")
ranef(fit.mm.nest.fwREML)
anova(fit.mm.nest.fwREML) # This test does also not take random effects into account.


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
fit.mm.autocor <- lme(clo ~ tOut,
                      random = ~1|subDay,
                      correlation = corAR1(form = ~time|subDay),
                        # corAR1(form = ~1|subDay/day), 
                      data = c.data, method = "ML")

add1(object = fit.mm.autocor, scope = ~.+ tInOp + time + sex, test = "Chisq")
# add sex
fit.mm.autocor1 <- update(fit.mm.autocor, .~.+ sex+tInOp)
anova(fit.mm.autocor1, fit.mm.autocor) # new model is better

add1(object = fit.mm.autocor1, scope = ~.+  time + tOut*tInOp*time*sex, test = "Chisq")
# add tInOp*sex
fit.mm.autocor2 <- update(fit.mm.autocor1, .~.+ tOut*sex + sex*tInOp)
anova(fit.mm.autocor2, fit.mm.autocor1) # new model is better

add1(object = fit.mm.autocor2, scope = ~.+ tInOp + time + tOut*tInOp*time*sex + I(tOut^2) + I(tInOp^2) + I(time^2), test = "Chisq")
# Nothing more to add


# final model:
fit.mm.autocor2$terms

# The forward model is better and simpler
fit.mm.autocor <- lme(clo ~ tOut,
                      random = ~1|subDay,
                      correlation = corAR1(form = ~time|subDay),
                      data = c.data, method = "REML")
ranef(fit.mm.autocor2.fwREML)
anova(fit.mm.autocor2.fwREML) # This test does also not take random effects into account.


# Final models for later comparison
final.model.REML3 <- fit.mm.autocor2.fwREML
final.model.ML3 <- fit.mm.autocor2



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

par(mfrow = c(2, 1))
acf(residuals(final.model.ML3, retype="normalized"), main = "Auto-correlation for final.model.ML3 (ACF)")
pacf(residuals(final.model.ML3, retype="normalized"), main = "Auto-correlation for final.model.ML3 (PACF)")


plot(Variogram(final.model.ML3))
