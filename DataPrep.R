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

#Remove rows containing NA
dioxin <- drop_na(dioxin, names(dioxin))
dioxin$logDiox <- log(dioxin$DIOX)

#Block effects: PLANT (3 plants, RENO_N, RENO_S and KARA), TIME (For RENO_N the experiment
#was repeated at a later time point, 2, as well.), LAB (Two labs. One in DK and one in USE)
#considerable measurement noise is expected.
dioxin$LOAD_Ordinal <- rep(3, dim(dioxin)[1])
dioxin$LOAD_Ordinal     <- dioxin$LOAD_Ordinal - 2*as.numeric(dioxin$LOAD == "L")
dioxin$LOAD_Ordinal     <- dioxin$LOAD_Ordinal - as.numeric(dioxin$LOAD == "N")

dioxin$OXYGEN_Ordinal <-  rep(3, dim(dioxin)[1])
dioxin$OXYGEN_Ordinal     <- dioxin$OXYGEN_Ordinal - 2*as.numeric(dioxin$OXYGEN == "L")
dioxin$OXYGEN_Ordinal     <- dioxin$OXYGEN_Ordinal -   as.numeric(dioxin$OXYGEN == "N")

dioxin$PRSEK_Ordinal <-  rep(3, dim(dioxin)[1])
dioxin$PRSEK_Ordinal     <- dioxin$PRSEK_Ordinal - 2*as.numeric(dioxin$PRSEK == "L")
dioxin$PRSEK_Ordinal     <- dioxin$PRSEK_Ordinal -   as.numeric(dioxin$PRSEK == "N")

dioxin$PLANT_RENO_N <- as.numeric(dioxin$PLANT == "RENO_N") #0 means Reno_S
dioxin$PLANT_KARA   <- as.numeric(dioxin$PLANT == "KARA")
dioxin$LAB_USA_or_KK      <- as.numeric(dioxin$LAB == "USA")