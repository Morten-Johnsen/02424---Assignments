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
  select(-PLANT, -LAB, -OXYGEN, -LOAD, -PRSEK, -OBSERV) %>%
  melt() %>%
  ggplot()+
  geom_histogram(aes(x = value), bins = 10)+
  facet_wrap(~variable, scales = "free")

#Block effects: PLANT (3 plants, RENO_N, RENO_S and KARA), TIME (For RENO_N the experiment
#was repeated at a later time point, 2, as well.), LAB (Two labs. One in DK and one in USE)
#considerable measurement noise is expected.

dioxin %>%
  select(logDiox, TIME,# PLANT_RENO_S - is 0 in PLANT_RENO_N
         , OXYGEN_Ordinal
         , LOAD_Ordinal
         , PRSEK_Ordinal
         , O2, O2COR, NEFFEKT, QRAT) %>%
  ggpairs()

# Plot ordinal values versus the actually measured values
dioxin %>%
  select(OXYGEN_Ordinal
         , LOAD_Ordinal
         , PRSEK_Ordinal
         , O2, O2COR, NEFFEKT, QRAT) %>%
  ggpairs()

# Plot the values used in the first models: 
# Block values and the active values (ordinal)
dioxin %>%
  select(logDiox, TIME, LAB_USA_or_KK
         , PLANT_RENO_N, PLANT_KARA # PLANT_RENO_S - is 0 in PLANT_RENO_N
         , OXYGEN_Ordinal
         , LOAD_Ordinal
         , PRSEK_Ordinal) %>%
  ggpairs()


# Plot the values used in the second models: 
# Block values and the active values (measured)
dioxin %>%
  select(logDiox, TIME, factor(LAB)
         , PLANT
         , O2, O2COR, NEFFEKT, QRAT) %>%
  ggpairs()

#Active variables: 
#   "Theoretical": OXYGEN, LOAD, PRSEK.
#   "Measured": O2 (O2COR), NEFFEKT, QRAT


# Model(s) ----------------------------------------------------------------
#### 2) ####
# Model with only active and the block variables.
fit0 <- lm(logDiox ~ factor(OXYGEN) + LOAD_Ordinal + PRSEK_Ordinal + 
             factor(PLANT) + TIME + factor(LAB), data = dioxin)
summary(fit0)

# Model with only active and the block variables.
# -PRSEK
fit1 <- lm(logDiox ~ factor(OXYGEN) + LOAD_Ordinal  + 
             factor(PLANT) + TIME + factor(LAB), data = dioxin)
summary(fit1)

anova(fit0, fit1)
# Not significant change. Therefore, we like fit1 better because it's simple.

# -PRSEK
fit2 <- lm(logDiox ~  LOAD_Ordinal  + 
             factor(PLANT) + TIME + factor(LAB), data = dioxin)
summary(fit2)

anova(fit1,fit2)
# Significant change, thus, we need to stick with fit1. 

# Fit plot with labelled outliers
par(mfrow=c(2,2))
plot(fit1, pch = 19, col = 'gray50', 
     main = "MODEL VARIATION", 
     sub='green = variation accounted for by the model')


#### 3) ####
# Model with measured and the block variables.
fit3_0 <- lm(logDiox ~ O2COR + NEFFEKT + QRAT + 
             factor(PLANT) + TIME + factor(LAB), data = dioxin)
summary(fit3_0)

# Reduce model
# - QRAT
fit3_1 <- lm(logDiox ~ O2COR + NEFFEKT +
               factor(PLANT) + TIME + factor(LAB), data = dioxin)
summary(fit3_1)

# Not significant change. Therefore, we like fit3_1 better because it's simple.
anova(fit3_0, fit3_1)

# Reduce model
# - PLANT
fit3_2 <- lm(logDiox ~ O2COR + NEFFEKT + TIME + factor(LAB), data = dioxin)
summary(fit3_2)

anova(fit3_1,fit3_2)
# Significant change! Therefore, we keep fit3_1

# Fit plot with labelled outliers
par(mfrow=c(2,2))
plot(fit3_1, pch = 19, col = 'gray50', 
     main = "MODEL VARIATION", 
     sub='green = variation accounted for by the model')

plot(fit1, pch = 19, col = 'gray50', 
     main = "MODEL VARIATION", 
     sub='green = variation accounted for by the model')

# We see that the qq plots have a better fit for the model with measured variables (fit2)


predict(fit3_1,newdata=data.frame(PLANT="RENO_N", O2COR = 0.5, NEFFEKT = -0.01, 
                                  TIME = factor("1"), LAB = "KK"),interval="prediction")


confint(fit3_1)


dioxin %>%
  select(logDiox, TIME,# PLANT_RENO_S - is 0 in PLANT_RENO_N
         , QROEG,TOVN,TROEG,POVN,CO2,logCO,SO2,logHCL,H2O) %>%
  ggpairs()


add1(fit3_1, scope=~.+QROEG+TOVN+TROEG+POVN+CO2+logCO+SO2+logHCL+H2O,test="F")

add1(fit3_1, scope=~.+QROEG+TOVN+TROEG+POVN+CO2+logCO+SO2+logHCL+H2O+I(QROEG^2)+I(TOVN^2)+
       I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(logCO^2)+I(SO2^2)+I(logHCL^2)+I(H2O^2)+
       I(QROEG*TOVN)+I(QROEG*TROEG)+I(QROEG*POVN)+I(QROEG*SO2)+I(QROEG*H2O)+
       I(TOVN*CO2)+I(TROEG*logCO)+I(TROEG*SO2)+I(TROEG*logHCL)+I(TROEG*H2O)+
       I(CO2*H2O)+I(SO2*H2O)
       ,test="F")

fit7 <- update(fit3_1, .~. + logHCL,data=dioxin)
summary(fit7)

add1(fit7, scope=~.+I(logHCL*TROEG)
     ,test="F")

fit7_1 <- update(fit7, .~. + I(logHCL*TROEG),data=dioxin)
summary(fit7_1)

anova(fit7,fit7_1)


fit7_2 <- update(fit7_1, .~. + CO2,data=dioxin)
summary(fit7_2)

anova(fit7_1,fit7_2)

fit7_3 <- update(fit7_2, .~. + POVN,data=dioxin)
summary(fit7_3)

anova(fit7_2,fit7_3)

plot(fit7_3)

# 9)
1/sd(dioxin[dioxin$LAB=="KK",]$logDiox)
1/sd(dioxin[dioxin$LAB=="USA",]$logDiox)
