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
earinfect$healthy <- earinfect$persons - earinfect$infections


# Plots -------------------------------------------------------------------

hist(earinfect$infections)
ggpairs(earinfect)

ggplot(earinfect) + 
  geom_histogram(aes(x = infections, y = after_stat(density)), bins = 15,  colour = "white", position = "identity", alpha = 0.4) +
  theme_bw() +
  labs(y = "", colour = "Distribution", title = "Earinfections") 
ggsave("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 2/earinfect.png", width = 20, height = 10, units = "cm")

mean(earinfect$infections)
var(earinfect$infections)


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




# Quasipoisson - for over-dispersion


# Plots -------------------------------------------------------------------

1-pchisq(fit.pois$deviance, df = fit.pois$df.residual, lower.tail = F) #model is appropriate
par(mfrow = c(2,2))
plot(fit.pois)

1-pchisq(fit.lognorm$deviance, df = fit.lognorm$df.residual, lower.tail = F) #model is appropriate
plot(fit.lognorm)

# Predict
plot(log(fitted(fit.pois)),log((earinfect$infections-fitted(fit.pois))^2),
     xlab=expression(hat(mu)),ylab=expression((y-hat(mu))^2),
     pch=20,col="blue", ylim = c(-4, 4))
abline(0,1) # mean
# We can see that the majority of the variance is smaller than the mean, 
# which is a warning of overdispersion?


# Dispersion parameter being larger than 1 - hence overdispersion
dp = sum(residuals(fit.pois,type ="pearson")^2)/fit.pois$df.residual
dp

# Dispersion test
library(AER)
dispersiontest(fit.pois)
# This overdispersion test reports the significance of the overdispersion issue within the model
# Eg. it is not significant


# Check how the parameters are affected of the disperiosn
summary(fit.pois,dispersion = dp)
# Not quite affected - a lot of parameters are still significant



# Extract the difference between the predicted and actual values
resDev <- residuals(fit.pois, type="deviance")
resDevPear <- residuals(fit.pois, type = "pearson")

# Residual analysis of the variation in cloting insulation level between females and men.
par(mfrow = c(1, 1))
plot(jitter(as.numeric(earinfect$sex), amount=0.1), resDev,
     xlab="Sex", ylab="Deviance residuals", cex=0.6,
     axes=FALSE)
box()
axis(1,label=c("Female", "Male" ),at=c(1,2))
axis(2)
# The deviance is quite distributed evenly for the sex, but we do see a larger variance for females

confint(fit.pois)

par(mfrow = c(1, 2))
plot(earinfect$infections, resDev, xlab='infections', ylab='Deviance residuals')
plot(earinfect$persons, resDev, xlab='persons', ylab='Deviance residuals')
# A bit random, which is good

#gender-specific residual analysis
earinfect$pred <- predict(fit.pois)
earinfect$pearson <- residuals(fit.pois, type = "pearson")

#Check whether residuals are i.i.d.
ggplot(earinfect)+
  geom_boxplot(aes(x = sex, y = pearson, fill = sex))+
  theme_bw()+
  ggtitle("Residual variation difference between genders")
#=> Residuals are ish identically distributed. This will be addressed later.


basicPlot <- function(...){
  plot(infections ~ age + sex + location + swimmer + age:sex + 
         age:location + sex:swimmer + location:swimmer + age:sex:location + 
         offset(persons), data=earinfect, bty="n", lwd=2,
       main="Number of earinfections", col="#00526D", 
       xlab="Infections", 
       ylab="?", ...)
  axis(side = 1, col="grey")
  axis(side = 2, col="grey")
}
temp <- 1:24
p.pois <- predict(fit.pois, type="response")
basicPlot()
points(temp, p.pois, col="orange", lwd=2)
legend(x="topleft", 
       legend=c("observation", 
                "Poisson (log) GLM"),
       col=c("#00526D","orange"),  
       bty="n", lwd=rep(2,5))


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
