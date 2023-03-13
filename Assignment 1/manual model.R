rm(list=ls())
library(GGally)
library(corrplot)
library(car)
library(MASS)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/")
}

source("DataPrep.R")

y <- dioxin$logDiox

dioxin$PLANT
dioxin$TIME
dioxin$LAB
#include lab = KK, PLANT = kara and time = 1 in the intercept estimate and model the other

RENO_N <- as.numeric(dioxin$PLANT == "RENO_N")
RENO_S <- as.numeric(dioxin$PLANT == "RENO_S")
TIME_2 <- as.numeric(dioxin$TIME == 2)
LAB_USA <- as.numeric(dioxin$LAB == "USA")

#create design matrix
X <- cbind(1, dioxin$O2COR, dioxin$NEFFEKT, RENO_N, RENO_S, TIME_2, LAB_USA)

beta <- solve(t(X)%*%X)%*%t(X)%*%y #p. 49 eq. 3.27
library(numDeriv)
sd_e <- sqrt(t(y - X %*% beta) %*% (y - X %*% beta) / (dim(X)[1] - dim(X)[2])) #p. 53 eqe. 3.40

get.hessian <- function(beta, sigma = sd_e, design = X, target = y){
  e <- target - design %*% beta
  return(-sum(dnorm(e, sd = sigma, log = T)))
}

sqrt(diag(solve(hessian(get.hessian, x = beta))))

