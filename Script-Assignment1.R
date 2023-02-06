rm(list=ls())
library(tidyverse)
library(reshape2)
library(gridExtra)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/")
}

dioxin <- read.csv("./dioxin.csv")

tibble(dioxin)
