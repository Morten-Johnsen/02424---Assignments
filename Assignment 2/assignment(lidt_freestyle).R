## Assignment 2 - part A
library(ggplot2)
library(tidyverse)

# Load data ---------------------------------------------------------------
if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 2/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 2")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 2/")
}

data <- read.table("clothing.csv", 
                   header = TRUE, sep = ",")
str(data)
data <- tibble(data)
data$SEX <- ifelse(data$sex == "female", 1, 0) # Convert sex to binary
data$resp <- cbind(data$clo*100, rep(100, length(data$clo))-data$clo*100)


# Find sutiable GLM for clothing insulation level and compare the results
# Should only concern the modelling of clothing insulation level based on
# indoor and outdoor temperature
data1 <- data %>%
  select(-subjId, -day, -time)
data1$sex <- as.factor(data1$sex)

# 'clo' clothing level is between 0 and 1
# The other variables for temperature are continuos (not binomial)

# Distribution of sex
table(data1$sex)


# Plots -------------------------------------------------------------------
par(mfrow = c(1, 3))
plot(clo ~ tOut, col="green", data=data1)
plot(clo ~ tInOp, col="green", data=data1)
plot(clo ~ sex, col="green", data=data1)
# 
# p = data$clo
# logit <- log(p/(1-p))
# plot(data$clo, logit, las=1,
#      xlab='clo', ylab="Logit(clo)")

# Models -------------------------------------------------------------------

# All variables combined
fit0 <- glm(formula = resp ~ tOut*tInOp*sex,
            data = data1,
            family = binomial(link="logit"))
summary(fit0) # None is significant
drop1(fit0, test = "Chisq")

fit1 <- update(fit0,.~.-tOut:tInOp:sex)
summary(fit1)
drop1(fit1, test = "Chisq")

fit2 <- update(fit1,.~.-tOut:tInOp-tOut:sex-tInOp:sex)
summary(fit2) # Somthing are significant now
drop1(fit2, test = "Chisq")

# fit2 is the same as a very simple additive model
fit3 <- update(fit2,.~.-tInOp)
summary(fit3) # Somthing are significant now
drop1(fit3, test = "Chisq")

anova(fit0,fit1,test="Chisq")
anova(fit0,fit2,test="Chisq")
anova(fit0,fit3,test="Chisq")

# We see that none of the anova tests are significant, meaning that we will take/Accept the simpler model(s)
fit3 = fit0

# Residuals ---------------------------------------------------------------

# Extract the difference between the predicted and actual values
resDev <- residuals(fit3, type="deviance")
resDevPear <- residuals(fit3, type = "pearson")

# Residual analysis of the variation in cloting insulation level between females and men.
par(mfrow = c(1, 1))
plot(jitter(as.numeric(as.factor(data$sex)), amount=0.1), resDev,
     xlab="Sex", ylab="Deviance residuals", cex=0.6,
     axes=FALSE)
box()
axis(1,label=c("Female", "Male" ),at=c(1,2))
axis(2)
# The deviance is quite distributed evenly for the sex, but we do see a larger variance for females

confint(fit3)

par(mfrow = c(1, 3))
plot(data1$clo, resDev, xlab='Clothing insulation level', ylab='Deviance residuals')
plot(data1$tOut, resDev, xlab='Temperature outdoot', ylab='Deviance residuals')
plot(data1$tInOp, resDev, xlab='Operating temperature indoor', ylab='Deviance residuals')
# We do see that the residuals for 'clo' is linear distributed. We want more randomness...

par(mfrow = c(2, 2))
plot(fit3)

# Cooks distance
cooksD  <- cooks.distance(fit3)
cooksD[cooksD>(3*mean(cooksD))]
influential <- data1[cooksD>(3*mean(cooksD)),]



# New model with 'subjId' --------------------------------------------------------------

data2 <- data %>%
  select(-sex, -day, -time)
data2$subjId <- as.factor(data2$subjId)

par(mfrow = c(1, 1))
plot(clo ~ subjId, col="green", data=data2)

# Fit
# All variables combined
fit0 <- glm(formula = resp ~ tOut*tInOp*subjId, #clo
            data = data2,
            # family = gaussian(link = inverse))
            family = binomial(link="logit"))
summary(fit0) # None is significant
drop1(fit0, test = "Chisq")

fit1 <- update(fit0,.~.-tOut:tInOp:subjId)
summary(fit1)
drop1(fit1, test = "Chisq")

fit2 <- update(fit1,.~.-tOut:tInOp-tOut:subjId-tInOp:subjId)
summary(fit2) # Somthing are significant now
drop1(fit2, test = "Chisq")

# fit2 is the same as a very simple additive model
fit3 <- update(fit2,.~.-subjId)#-tInOp)
summary(fit3) # Somthing are significant now
drop1(fit3, test = "Chisq")

anova(fit0,fit1,test="Chisq")
anova(fit0,fit2,test="Chisq")
anova(fit0,fit3,test="Chisq")

