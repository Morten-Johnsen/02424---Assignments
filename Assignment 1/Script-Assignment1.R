rm(list=ls())
library(GGally)
library(corrplot)
library(car)
library(MASS)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 1/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 1")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 1/")
}

source("DataPrep.R")

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

# Histogram of boxcox and log-dioxin
dioxin$DIOX_boxcox <- Trans.eq1(lambda.hat, dioxin$DIOX)
hist(dioxin$DIOX_boxcox, breaks = 5)

# Histogram of variables
postscript("histograms_for_all.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 10, width = 10)
dioxin %>%
  dplyr::select(-PLANT, -LAB, -OXYGEN, -LOAD, -PRSEK, -OBSERV, -DIOX_boxcox, -logCO) %>%
  melt() %>%
  # Color code passive vs active, and measured vs. ordinal
  mutate(color = c(rep("black", 52), 
                   rep("green", 52*4), 
                   rep("orange", 52*9), 
                   rep("black", 52), 
                   rep("orange", 52*2), 
                   rep("red", 52*3))) %>%
  ggplot()+
  geom_histogram(aes(x = value, fill = color), bins = 10)+
  facet_wrap(~variable, scales = "free")+
  theme(legend.position = "none")
dev.off()

# orange = passive (QROEG, TOVN, TROEG, POVN, CO2, CO, SO2, HCL, H2O)
# blue = active    (OXYGEN, LOAD, PRSEK)
# green = measured (O2, O2COR, NEFFEKT, QRAT)
# red = ordinal    (PALNT, TIME, LAB)
# black = DIOX


# Take out tables of data that can't be plotted
table(dioxin$PLANT)
table(dioxin$LAB)
table(dioxin$TIME)


# Make correlation plot of all variables:
#       - Ordinal vs. measured
#       - Variables in first model
#       - Variables in second model
par(mfrow=c(1,1))
postscript("corplot_ordinal_vs_measured.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 10, width = 10)
dioxin %>%
  dplyr::select(logDiox,
                O2, O2COR, NEFFEKT, QRAT
                #, TIME#, LAB
                #, PLANT # PLANT_RENO_S - is 0 in PLANT_RENO_N
                , OXYGEN_Ordinal
                , LOAD_Ordinal
                , PRSEK_Ordinal) %>%
  cor() %>%
  corrplot(method = 'color')
# corrplot()
dev.off()

# First model
postscript("ggpairs_first_model.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 10, width = 10)
dioxin %>%
  dplyr::select(logDiox, # PLANT_RENO_S - is 0 in PLANT_RENO_N
                PLANT,
                TIME,
                LAB,
                LOAD_Ordinal,
                OXYGEN_Ordinal) %>%
  ggpairs()
dev.off()


# Second model
postscript("ggpairs_second_model.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 10, width = 10)
dioxin %>%
  dplyr::select(logDiox, # PLANT_RENO_S - is 0 in PLANT_RENO_N
                PLANT,
                TIME,
                LAB,
                NEFFEKT,
                O2COR) %>%
  ggpairs()
dev.off()


#Active variables: 
#   "Theoretical": OXYGEN, LOAD, PRSEK.
#   "Measured": O2 (O2COR), NEFFEKT, QRAT

# Model(s) ----------------------------------------------------------------
names(dioxin)

#Block effects: PLANT (3 plants, RENO_N, RENO_S and KARA), TIME (For RENO_N the experiment
#was repeated at a later time point, 2, as well.), LAB (Two labs. One in DK and one in USA)
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

confint(fit2)

# Check outliers
cooksD  <- cooks.distance(fit2)
cooksD[cooksD>(3*mean(cooksD))]
influential <- dioxin[cooksD>(3*mean(cooksD)),]
mean(dioxin$DIOX)

postscript("residualplotsFirstModel.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 10, width = 11)
par(mfrow=c(2,2))
plot(fit2,cex.lab=1.3, cex.axis=1.3, cex.main=2, cex.sub=2)
dev.off()

# Make figures that highlight the outliers
# First with diox
postscript("outliers_diox.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 4, width = 13)
plotData <- dioxin %>%
  dplyr::select(DIOX, PLANT, TIME, LAB, LOAD, OXYGEN)  %>%
  melt(id = 'DIOX') 
ggplot(plotData, aes(x = value, y = DIOX)) +
  geom_point() +
  xlab("") + 
  geom_point(data = plotData[cooksD>(3*mean(cooksD)),], size = 3, color = "red") +
  facet_grid(. ~ variable,scales = "free_x", switch = 'x') +
  theme(text = element_text(size = 18))    
dev.off()

# Then with logdiox
postscript("outliers_logdiox.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 4, width = 13)
plotData1 <- dioxin %>%
  dplyr::select(logDiox, PLANT, TIME, LAB, LOAD, OXYGEN)  %>%
  melt(id = 'logDiox') 
ggplot(plotData1, aes(x = value, y = logDiox)) +
  geom_point() + 
  xlab("") + 
  geom_point(data = plotData1[cooksD>(3*mean(cooksD)),], size = 3, color = "red") +
  facet_grid(. ~ variable,scales = "free_x", switch = 'x') + 
  theme(text = element_text(size = 18))    
dev.off()

#### 3) ####
fit_obs <- lm(logDiox ~ O2COR + NEFFEKT + QRAT + PLANT + TIME + LAB, data = dioxin)
summary(fit_obs)
Anova(fit_obs, type = "III")

drop1(fit_obs, test = "F")
fit_obs1 <- update(fit_obs, .~.-QRAT)
drop1(fit_obs1, test = "F")


postscript("residualplotsSecondModel.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 10, width = 11)
par(mfrow=c(2,2))
plot(fit_obs1,cex.lab=1.3, cex.axis=1.3, cex.main=2, cex.sub=2)
dev.off()

anova(fit_obs, fit_obs1) #model performance are the same
Anova(fit_obs1, type = "III") #type 3 anova

summary(fit_obs1)

# Check outliers
cooksD  <- cooks.distance(fit_obs1)
cooksD[cooksD>(3*mean(cooksD))]
influential <- dioxin[cooksD>(3*mean(cooksD)),]
mean(dioxin$DIOX)

# Make figures that highlight the outliers

# First with diox
postscript("outliers_diox_measured.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 4, width = 13)
plotData <- dioxin %>%
  dplyr::select(DIOX, PLANT, TIME, LAB, LOAD, OXYGEN)  %>%
  melt(id = 'DIOX') 
ggplot(plotData, aes(x = value, y = DIOX)) +
  geom_point() +
  xlab("") + 
  geom_point(data = plotData[cooksD>(3*mean(cooksD)),], size = 3, color = "red") +
  facet_grid(. ~ variable,scales = "free_x", switch = 'x') +
  theme(text = element_text(size = 18))    
dev.off()

# Then with logdiox
postscript("outliers_logdiox_measured.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 4, width = 13)
plotData1 <- dioxin %>%
  dplyr::select(logDiox, PLANT, TIME, LAB, LOAD, OXYGEN)  %>%
  melt(id = 'logDiox') 
ggplot(plotData1, aes(x = value, y = logDiox)) +
  geom_point() + 
  xlab("") + 
  geom_point(data = plotData1[cooksD>(3*mean(cooksD)),], size = 3, color = "red") +
  facet_grid(. ~ variable,scales = "free_x", switch = 'x') + 
  theme(text = element_text(size = 18))    
dev.off()


#### 4) ####
# Make prediction
predict(fit_obs1, newdata = data.frame("PLANT" = factor("RENO_N"), 
                                       "LAB" = factor("KK"),
                                       "TIME" = factor(1),
                                       "O2COR" = 0.5,
                                       "NEFFEKT" = -0.01), 
        interval = "prediction",
        level = 0.95)

#### 5) ####
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
# Look at the passive variables to choose some interactions to model
dioxin %>%
  dplyr::select(logDiox # PLANT_RENO_S - is 0 in PLANT_RENO_N
         , QROEG,TOVN,TROEG,POVN,CO2,boxcoxCO,SO2,logHCL,H2O) %>%
  ggpairs()


# Look at different additions to the model
add1(fit_obs1, scope=~.+QROEG+TOVN+TROEG+POVN+CO2+boxcoxCO+SO2+logHCL+H2O+I(QROEG^2)+I(TOVN^2)+
       I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(boxcoxCO^2)+I(SO2^2)+I(logHCL^2)+I(H2O^2)+
       I(QROEG*TOVN)+I(QROEG*TROEG)+I(QROEG*POVN)+I(QROEG*SO2)+I(QROEG*H2O)+
       I(TOVN*CO2)+I(TROEG*boxcoxCO)+I(TROEG*SO2)+I(TROEG*logHCL)+I(TROEG*H2O)+
       I(CO2*H2O)+I(SO2*H2O)
     ,test="F")

# Add logHCL to model, because it is the most significant addition
fit_pas1 <- update(fit_obs1, .~. + logHCL,data=dioxin)
summary(fit_pas1)

anova(fit_obs1,fit_pas1)

add1(fit_pas1, scope=~.+QROEG+TOVN+TROEG+POVN+CO2+boxcoxCO+SO2+logHCL+H2O+I(QROEG^2)+I(TOVN^2)+
       I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(boxcoxCO^2)+I(SO2^2)+I(logHCL^2)+I(H2O^2)+
       I(QROEG*TOVN)+I(QROEG*TROEG)+I(QROEG*POVN)+I(QROEG*SO2)+I(QROEG*H2O)+
       I(TOVN*CO2)+I(TROEG*boxcoxCO)+I(TROEG*SO2)+I(TROEG*logHCL)+I(TROEG*H2O)+
       I(CO2*H2O)+I(SO2*H2O)
     ,test="F")

# Add TROEG to model, because it is the most significant addition
fit_pas2 <- update(fit_pas1, .~. + TROEG,data=dioxin)
summary(fit_pas2)
# Check difference in models
anova(fit_pas1,fit_pas2)
# Models are significantly different. Therefore, we keep the new one.

add1(fit_pas2, scope=~.+QROEG+TOVN+TROEG+POVN+CO2+boxcoxCO+SO2+logHCL+H2O+I(QROEG^2)+I(TOVN^2)+
       I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(boxcoxCO^2)+I(SO2^2)+I(logHCL^2)+I(H2O^2)+
       I(QROEG*TOVN)+I(QROEG*TROEG)+I(QROEG*POVN)+I(QROEG*SO2)+I(QROEG*H2O)+
       I(TOVN*CO2)+I(TROEG*boxcoxCO)+I(TROEG*SO2)+I(TROEG*logHCL)+I(TROEG*H2O)+
       I(CO2*H2O)+I(SO2*H2O)
     ,test="F")

# Add CO2 to model, because it is the most significant addition
fit_pas3 <- update(fit_pas2, .~. + CO2,data=dioxin)
summary(fit_pas3)
# Check difference in models
anova(fit_pas2,fit_pas3)
# Models are significantly different. Therefore, we keep the new one.

add1(fit_pas3, scope=~.+QROEG+TOVN+TROEG+POVN+CO2+boxcoxCO+SO2+logHCL+H2O+I(QROEG^2)+I(TOVN^2)+
       I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(boxcoxCO^2)+I(SO2^2)+I(logHCL^2)+I(H2O^2)+
       I(QROEG*TOVN)+I(QROEG*TROEG)+I(QROEG*POVN)+I(QROEG*SO2)+I(QROEG*H2O)+
       I(TOVN*CO2)+I(TROEG*boxcoxCO)+I(TROEG*SO2)+I(TROEG*logHCL)+I(TROEG*H2O)+
       I(CO2*H2O)+I(SO2*H2O)
     ,test="F")

# Add POVN to model, because it is the most significant addition
fit_pas4 <- update(fit_pas3, .~. + POVN,data=dioxin)
summary(fit_pas4)
# Check difference in models
anova(fit_pas3,fit_pas4)
# Models are significantly different. Therefore, we keep the new one.

add1(fit_pas4, scope=~.+QROEG+TOVN+TROEG+POVN+CO2+boxcoxCO+SO2+logHCL+H2O+I(QROEG^2)+I(TOVN^2)+
       I(TROEG^2)+I(POVN^2)+I(CO2^2)+I(boxcoxCO^2)+I(SO2^2)+I(logHCL^2)+I(H2O^2)+
       I(QROEG*TOVN)+I(QROEG*TROEG)+I(QROEG*POVN)+I(QROEG*SO2)+I(QROEG*H2O)+
       I(TOVN*CO2)+I(TROEG*boxcoxCO)+I(TROEG*SO2)+I(TROEG*logHCL)+I(TROEG*H2O)+
       I(CO2*H2O)+I(SO2*H2O)
     ,test="F")

# Do not add any more terms.

postscript("residualplotsFinalModel.eps", horizontal = FALSE, onefile = FALSE, paper = "special",height = 10, width = 11)
par(mfrow=c(2,2))
plot(fit_pas4,cex.lab=1.3, cex.axis=1.3, cex.main=2, cex.sub=2)
dev.off()
plot(fit_pas4)
confint(fit_pas4)

# Leverages
cooksD  <- cooks.distance(fit_pas4)
cooksD[cooksD>(3*mean(cooksD))]
influential <- dioxin[cooksD>(3*mean(cooksD)),]



# possible likelihood for 9)
ml.function <- function(theta){
  return(-sum(dnorm(filter(dioxin, LAB == "KK")$logDiox, mean = theta[1], sd = theta[2], log = T)))
}

est.KK <- nlminb(c(0,1), objective = ml.function)$par
est.KK

ml.function2 <- function(theta){
  return(-sum(dnorm(filter(dioxin, LAB == "USA")$logDiox, mean = theta[1], sd = theta[2], log = T)))
}
est.USA <- nlminb(c(0,1), objective = ml.function2)$par
est.USA

1/est.KK[2]
1/est.USA[2]

y <- dioxin$logDiox
RENO_N <- as.numeric(dioxin$PLANT == "RENO_N")
RENO_S <- as.numeric(dioxin$PLANT == "RENO_S")
TIME_2 <- as.numeric(dioxin$TIME == 2)
LAB_USA <- as.numeric(dioxin$LAB == "USA")

#create design matrix
X <- cbind(1, dioxin$O2COR, dioxin$NEFFEKT, RENO_N, RENO_S, TIME_2, LAB_USA, 
           dioxin$logHCL, dioxin$CO2)

# Make function that re-estimate beta-parameters together with the weights.
ml.function <- function(theta, y, X){
  wUSA <- theta[1]
  wKK <- theta[2]
  w.diag <- rep(NA, dim(X)[1])
  w.diag[X[,7] == 1] <- wUSA
  w.diag[X[,7] == 0] <- wKK
  SIGMA <- diag(w.diag)
  
  n <- length(theta)
  beta <- matrix(theta[3:n], nrow = 9, ncol = 1)
  
  sd_e <<- as.numeric(sqrt(t(y - X %*% beta) %*% solve(SIGMA) %*% (y - X %*% beta) / (dim(X)[1] - dim(X)[2]))) #p. 53 eqe. 3.40
  return(-sum(log(1/(sd_e*w.diag * sqrt(2*pi)) * exp(-1/2 * ((y-X%*%beta)/(sd_e*w.diag))^2))))
}

theta.hat <- nlminb(start = c(1,1, rep(1,9)), objective = ml.function, y = y, X = X)
theta.hat$par
ml.function(theta.hat$par, y, X)
sd_e

uncertainties <- sqrt(diag(solve(hessian(func = ml.function, x = theta.hat$par, y = y, X = X))))
uncertainties

theta.hat$par/uncertainties
t(t(2*pt(abs(theta.hat$par/uncertainties), df = 52 - 11, lower.tail = F)))


w.diag <- rep(NA, dim(dioxin)[1])
w.diag[X[,7] == 0] <- theta.hat$par[2]
w.diag[X[,7] == 1] <- theta.hat$par[1]


fit <- lm(logDiox ~ O2COR + NEFFEKT + PLANT + TIME + LAB + logHCL + CO2 + POVN +TROEG , data=dioxin, weights = 1/w.diag^2)
logLik(fit)
coefficients(fit)
