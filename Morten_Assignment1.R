rm(list=ls())
library(GGally)
library(corrplot)
library(MASS)
library(car)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/")
}

source("DataPrep.R")

names(dioxin)

#Block effects: PLANT (3 plants, RENO_N, RENO_S and KARA), TIME (For RENO_N the experiment
#was repeated at a later time point, 2, as well.), LAB (Two labs. One in DK and one in USE)
#considerable measurement noise is expected.
#Active variables: 
#   "Theoretical": OXYGEN, LOAD, PRSEK.
#   "Measured": O2 (O2COR), NEFFEKT, QRAT
dioxin$plant <- factor(dioxin$PLANT)
dioxin$time <- factor(dioxin$TIME)
dioxin$lab <- factor(dioxin$LAB_USA_or_KK)
dioxin$oxygen <- factor(dioxin$OXYGEN)
dioxin$load <- factor(dioxin$LOAD)
dioxin$prsek <- factor(dioxin$PRSEK)


#### 2) ####
fit <- lm(logDiox ~ LOAD_Ordinal*OXYGEN_Ordinal*PRSEK_Ordinal + plant + time + lab, data = dioxin)
Anova(fit, type = "III")
summary(fit)
#plot(fit)
drop1(fit, test = "F")
fit2 <- update(fit, .~.-LOAD_Ordinal:OXYGEN_Ordinal:PRSEK_Ordinal)
drop1(fit2, test = "F")
fit3 <- update(fit2, .~.-LOAD_Ordinal:OXYGEN_Ordinal)
drop1(fit3, test = "F")
fit4 <- update(fit3, .~.-LOAD_Ordinal:PRSEK_Ordinal)
drop1(fit4, test = "F")
fit5 <- update(fit4, .~.-OXYGEN_Ordinal:PRSEK_Ordinal)
drop1(fit5, test = "F")
fit6 <- update(fit5, .~.-PRSEK_Ordinal)
drop1(fit6, test = "F")

anova(fit, fit6)
Anova(fit6, type = "III")

par(mfrow=c(2,2))
plot(fit6)

#### 3) ####
dioxin %>%
  dplyr::select(O2COR, NEFFEKT, QRAT, plant, time, lab, logDiox) %>%
  ggpairs()
fit_obs <- lm(logDiox ~ O2COR*NEFFEKT*QRAT + plant + time + lab, data = dioxin)
Anova(fit_obs, type = "III")

drop1(fit_obs, test = "F")
fit_obs1 <- update(fit_obs, .~.-O2COR:NEFFEKT:QRAT)
drop1(fit_obs1, test = "F")
fit_obs2 <- update(fit_obs1, .~.-O2COR:NEFFEKT)
drop1(fit_obs2, test = "F")
fit_obs3 <- update(fit_obs2, .~.-O2COR:QRAT)
drop1(fit_obs3, test = "F")
fit_obs4 <- update(fit_obs3, .~.-NEFFEKT:QRAT)
drop1(fit_obs4, test = "F")
fit_obs5 <- update(fit_obs4, .~.-QRAT)
drop1(fit_obs5, test = "F")

anova(fit_obs, fit_obs5) #model performance are the same
Anova(fit_obs5, type = "III") #type 3 anova

summary(fit_obs5)
plot(fit_obs5)
par(mfrow=c(1,1))
qqPlot(fit_obs5)
confint(fit_obs5, level = 0.99)

dioxin[c(13,24,11,20),] %>%
  dplyr::select(plant, time, lab, NEFFEKT, O2COR, logDiox)

dioxin$outlier <- "No"
dioxin$outlier[c(13,24,11,20)] <- "Yes"
dioxin$chrLog <- as.character(round(dioxin$logDiox,2))
library(gghighlight)
ggplot(dioxin)+
  geom_point(aes(x = OBSERV, y = logDiox, fill = outlier))+
  gghighlight(outlier == "Yes", label_key = chrLog)
