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
#### 2) ####
fit <- lm(logDiox ~ PLANT + TIME + LAB+LOAD_Ordinal + OXYGEN_Ordinal + PRSEK_Ordinal, data = dioxin)
#plot(fit)
drop1(fit, test = "F")
fit2 <- update(fit, .~.-PRSEK_Ordinal)
drop1(fit2, test = "F")

anova(fit, fit2)
Anova(fit2, type = "III")
summary(fit2)

par(mfrow=c(2,2))
plot(fit2)

#### 3) ####
fit_obs <- lm(logDiox ~ O2COR + NEFFEKT + QRAT + PLANT + TIME + LAB, data = dioxin)
Anova(fit_obs, type = "III")

drop1(fit_obs, test = "F")
fit_obs1 <- update(fit_obs, .~.-QRAT)
drop1(fit_obs1, test = "F")

anova(fit_obs, fit_obs1) #model performance are the same
Anova(fit_obs1, type = "III") #type 3 anova

summary(fit_obs1)
plot(fit_obs1)
par(mfrow=c(1,1))
qqPlot(fit_obs1)
confint(fit_obs1, level = 0.99)

dioxin[c(13,24,11,20),] %>%
  dplyr::select(PLANT, LAB, TIME, NEFFEKT, O2COR, logDiox)

dioxin$outlier <- "No"
dioxin$outlier[c(13,24,11,20)] <- "Yes"
dioxin$chrLog <- as.character(round(dioxin$logDiox,2))
library(gghighlight)
ggplot(dioxin)+
  geom_point(aes(x = PRSEK, y = logDiox, fill = outlier))+
  gghighlight(outlier == "Yes", label_key = chrLog)

#4)
predict(fit_obs1, newdata = data.frame("PLANT" = factor("RENO_N"), 
                                       "LAB" = factor("KK"),
                                       "TIME" = factor(1),
                                       "O2COR" = 0.5,
                                       "NEFFEKT" = -0.01), 
        interval = "prediction",
        level = 0.95)

#5)
AIC(fit2)
AIC(fit_obs1)
#Operating conditions make a difference.
coef(fit2)
coef(fit_obs1)
#dioxin emission depend on the operating conditions of oxygen (O2COR) and load (NEFFEKT)
#but not on PRSEK (QRAT). 
#Dioxin emission can be reduced by reducing oxygen and load.


#6)
#Differences between plants:
#PLANT_S predicts -2.11 log(ppm) less than plant KARA
#PLANT_N predicts -0.57 log(ppm) less than plant KARA

#7)
par(mfrow=c(2,3))
plot(dioxin$O2COR, fit_obs1$residuals)
plot(dioxin$NEFFEKT, fit_obs1$residuals)
plot(dioxin$TIME, fit_obs1$residuals)
plot(dioxin$PLANT, fit_obs1$residuals)
plot(dioxin$LAB, fit_obs1$residuals)


add1(fit_obs1, scope = ~.+O2COR*NEFFEKT+I(O2COR^2)+I(NEFFEKT^2)+QROEG+TOVN+TROEG+POVN+CO2+CO+SO2+HCL+O2COR+H2O+I(QROEG^2)+I(TOVN^2)+I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(CO^2)+I(SO2^2)+I(HCL^2)+I(H2O^2)
     , test = "F", data = dioxin)
model <- update(fit_obs1, .~.+HCL, data = dioxin)
add1(model, scope = ~.+I(O2COR^2)+I(NEFFEKT^2)+QROEG+TOVN+TROEG+POVN+CO2+CO+SO2+H2O+I(QROEG^2)+I(TOVN^2)+I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(CO^2)+I(SO2^2)+I(HCL^2)+I(H2O^2), test = "F", data = dioxin)
drop1(model, test = "F")

par(mfrow = c(2,2))
plot(model)


par(mfrow=c(1,1))
hist(dioxin$DIOX, breaks = 20)
