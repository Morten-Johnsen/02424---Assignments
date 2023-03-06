rm(list=ls())
library(tidyverse)
library(reshape2)
library(gridExtra)
library(numDeriv)
library(latex2exp)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/")
}

dioxin <- read.csv("./dioxin.csv")
tibble(dioxin)


# Function
Trans.eq1 <- function(lambda, y = dioxin$DIOX){
  y_lambda <- ((y)^lambda - 1)/lambda#, lambda > 0
  return(y_lambda)
}

#Remove rows containing NA
dioxin <- drop_na(dioxin, names(dioxin))
# Log transform
dioxin$logDiox <- log(dioxin$DIOX)
dioxin$logCO <- log(dioxin$CO)
dioxin$boxcoxCO <- Trans.eq1(-0.241, dioxin$logCO) #found through the boxcox script
dioxin$logHCL <- log(dioxin$HCL)

#Block effects: PLANT (3 plants, RENO_N, RENO_S and KARA), TIME (For RENO_N the experiment
#was repeated at a later time point, 2, as well.), LAB (Two labs. One in DK and one in USE)
#considerable measurement noise is expected.
dioxin$LOAD_Ordinal <- rep(1, dim(dioxin)[1])
dioxin$LOAD_Ordinal     <- dioxin$LOAD_Ordinal - 2*as.numeric(dioxin$LOAD == "L")
dioxin$LOAD_Ordinal     <- dioxin$LOAD_Ordinal - 1*as.numeric(dioxin$LOAD == "N")

dioxin$OXYGEN_Ordinal <-  rep(1, dim(dioxin)[1])
dioxin$OXYGEN_Ordinal     <- dioxin$OXYGEN_Ordinal - 2*as.numeric(dioxin$OXYGEN == "L")
dioxin$OXYGEN_Ordinal     <- dioxin$OXYGEN_Ordinal - 1*as.numeric(dioxin$OXYGEN == "N")

dioxin$PRSEK_Ordinal <-  rep(1, dim(dioxin)[1])
dioxin$PRSEK_Ordinal     <- dioxin$PRSEK_Ordinal - 2*as.numeric(dioxin$PRSEK == "L")
dioxin$PRSEK_Ordinal     <- dioxin$PRSEK_Ordinal - 1*as.numeric(dioxin$PRSEK == "N")

dioxin$PLANT <- factor(dioxin$PLANT)
dioxin$TIME <- factor(dioxin$TIME)
dioxin$LAB <- factor(dioxin$LAB)
dioxin$OXYGEN <- factor(dioxin$OXYGEN)
dioxin$LOAD <- factor(dioxin$LOAD)
dioxin$PRSEK <- factor(dioxin$PRSEK)






