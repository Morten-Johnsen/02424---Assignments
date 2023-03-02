rm(list=ls())
library(GGally)
library(corrplot)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/")
}

source("DataPrep.R")

Trans.eq1 <- function(lambda, y = dioxin$DIOX){
  y_lambda <- ((y)^lambda - 1)/lambda#, lambda > 0
  return(y_lambda)
}
lambda_NLL <- function(theta, y = dioxin$DIOX){
  lambda <- theta[1]
  y_lambda <- Trans.eq1(lambda, y)
  mu <- theta[2]#mean(y_lambda)
  sigma <- theta[3]#sd(y_lambda)
  NLL <- -sum(-1/2*log(sigma^2) - (y_lambda-mu)^2 / (2*sigma^2) + (lambda - 1)*log(y))
  return(NLL)
}
theta.hat <- nlminb(start=c(1, 0, 1), objective = lambda_NLL)
lambda.hat <- round(theta.hat$par[1], 3)
mu.hat <- round(theta.hat$par[2], 3)
sd.hat <- round(theta.hat$par[3], 3)

#See if the boxcox transformation is significantly different from the log-transformation
pf_lambda <- function(lambda, y = dioxin$DIOX){
  n <- length(y)
  y_lambda <- Trans.eq1(lambda)
  sd <- sd(y_lambda)
  ll <- -n/2 * log(sd^2) - n/2 + (lambda - 1)* sum(log(y))
  return(-ll)
}
lambda_interval <- seq(-.5,0.1,0.003)
pf_curve <- -sapply(lambda_interval, FUN = pf_lambda)-max(-sapply(lambda_interval, FUN = pf_lambda))
plot(lambda_interval
     ,pf_curve, type = "l", main = TeX("Profile Likelihood for \\lambda")
     ,ylim = c(-3,0))
grid()
abline(v = lambda.hat, lty = 2)
alpha <- 0.05
c <- -0.5 * qchisq(1-alpha, df = 1)
abline(h = c, col = "red")

#95% profilelikelihood CI for lambda include 0 and thus a log-transformation is sufficient and to be preferred over the boxcox.
#wald CI for confirmation:
sd_reg <- sqrt(diag(solve(hessian(func = lambda_NLL, x = theta.hat$par))))
lambda.hat + qt(c(alpha/2, 1-alpha/2), df = length(dioxin$DIOX) - 1) * sd_reg[1]


dioxin$DIOX_boxcox <- Trans.eq1(lambda.hat, dioxin$DIOX)
hist(dioxin$DIOX_boxcox, breaks = 5)

dioxin %>%
  dplyr::select(-PLANT, -LAB, -OXYGEN, -LOAD, -PRSEK, -OBSERV) %>%
  melt() %>%
  ggplot()+
  geom_histogram(aes(x = value), bins = 10)+
  facet_wrap(~variable, scales = "free")

#Block effects: PLANT (3 plants, RENO_N, RENO_S and KARA), TIME (For RENO_N the experiment
#was repeated at a later time point, 2, as well.), LAB (Two labs. One in DK and one in USE)
#considerable measurement noise is expected.

dioxin %>%
  dplyr::select(logDiox, TIME, LAB
         , PLANT # PLANT_RENO_S - is 0 in PLANT_RENO_N
         , OXYGEN_Ordinal
         , LOAD_Ordinal
         , PRSEK_Ordinal
         , O2, O2COR, NEFFEKT, QRAT) %>%
  ggpairs()

# Plot ordinal values versus the actually measured values
dioxin %>%
  dplyr::select(OXYGEN_Ordinal
         , LOAD_Ordinal
         , PRSEK_Ordinal
         , O2, O2COR, NEFFEKT, QRAT) %>%
  ggpairs()

# Plot the values used in the first models: 
# Block values and the active values (ordinal)
dioxin %>%
  dplyr::select(logDiox, TIME, LAB
         , PLANT # PLANT_RENO_S - is 0 in PLANT_RENO_N
         , OXYGEN_Ordinal
         , LOAD_Ordinal
         , PRSEK_Ordinal) %>%
  ggpairs()

dioxin %>%
  dplyr::select(logDiox#, TIME#, LAB
         #, PLANT # PLANT_RENO_S - is 0 in PLANT_RENO_N
         , OXYGEN_Ordinal
         , LOAD_Ordinal
         , PRSEK_Ordinal) %>%
  cor() %>%
  corrplot()

# Plot the values used in the second models: 
# Block values and the active values (measured)
dioxin %>%
  dplyr::select(logDiox, TIME, LAB
         , PLANT # PLANT_RENO_S - is 0 in PLANT_RENO_N
         , O2, O2COR, NEFFEKT, QRAT) %>%
  ggpairs()

dioxin %>%
  dplyr::select(logDiox#, TIME, LAB
         #, PLANT # PLANT_RENO_S - is 0 in PLANT_RENO_N
         , O2, O2COR, NEFFEKT, QRAT) %>%
  cor() %>%
  corrplot(method = 'ellipse', order = 'AOE', type = 'upper')


#Active variables: 
#   "Theoretical": OXYGEN, LOAD, PRSEK.
#   "Measured": O2 (O2COR), NEFFEKT, QRAT


# Model(s) ----------------------------------------------------------------
names(dioxin)

#Block effects: PLANT (3 plants, RENO_N, RENO_S and KARA), TIME (For RENO_N the experiment
#was repeated at a later time point, 2, as well.), LAB (Two labs. One in DK and one in USE)
#considerable measurement noise is expected.
#Active variables: 
#   "Theoretical": OXYGEN, LOAD, PRSEK.
#   "Measured": O2 (O2COR), NEFFEKT, QRAT
#### 2) ####
fit <- lm(logDiox ~ PLANT + TIME + LAB + LOAD_Ordinal + OXYGEN_Ordinal + PRSEK_Ordinal, data = dioxin)
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
#Differences between labs:
#LAB_USA register 0.4 log(ppm) higher than LAB_KK

ggplot(dioxin)+
  geom_boxplot(aes(x = PLANT, y =logDiox))

#7)
par(mfrow=c(2,3))
plot(dioxin$O2COR, fit_obs1$residuals)
plot(dioxin$NEFFEKT, fit_obs1$residuals)
plot(dioxin$TIME, fit_obs1$residuals)
plot(dioxin$PLANT, fit_obs1$residuals)
plot(dioxin$LAB, fit_obs1$residuals)

ggpairs(dioxin[,c(-(1:8),-26,-27)])

add1(fit_obs1, scope = ~.+O2COR*NEFFEKT+I(O2COR^2)+I(NEFFEKT^2)+QROEG+TOVN+TROEG+POVN+CO2+CO+SO2+HCL+O2COR+H2O+
       I(QROEG^2)+I(TOVN^2)+I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(CO^2)+I(SO2^2)+I(HCL^2)+I(H2O^2)+
       QROEG:TOVN+QROEG:POVN+QROEG:TROEG+logCO
     , test = "F", data = dioxin)
model <- update(fit_obs1, .~.+HCL, data = dioxin)
add1(model, scope = ~.+I(O2COR^2)+I(NEFFEKT^2)+QROEG+TOVN+TROEG+POVN+CO2+CO+SO2+H2O+I(QROEG^2)+I(TOVN^2)+I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(CO^2)+I(SO2^2)+I(HCL^2)+I(H2O^2), test = "F", data = dioxin)
drop1(model, test = "F")

par(mfrow = c(2,2))
plot(model)

par(mfrow=c(1,1))
hist(dioxin$DIOX, breaks = 20)
