rm(list=ls())

library(reshape)
library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 3/")
  figpath <- "/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 3/figs/"
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 3")
  figpath <- "~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 3/figs/"
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 3/")
  figpath <- "C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 3/figs"
}

clothing <- read.csv(file = "clothingFullAss03.csv", header = T)
head(clothing)
clothing$subjId <- factor(clothing$subjId)
#clothing$day <- factor(clothing$day)

unique(clothing$day)
unique(clothing$time)
c.data <- dplyr::select(clothing, -sex,-time2, -X, -subDay)
head(c.data)


melt(c.data,id = c('subjId','day') )%>%
  ggplot()+
  geom_point(aes(x = day, y=value,colour = subjId))+
  facet_wrap(~variable, scales = "free")+
  theme_bw()

melt(c.data,id = c('subjId','clo'))%>%
  ggplot()+
  geom_point(aes(x = value, y=clo,colour = subjId), show.legend = FALSE)+
  facet_wrap(~variable, scales = "free")+
  theme_bw()
