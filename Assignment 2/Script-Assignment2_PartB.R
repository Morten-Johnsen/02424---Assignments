rm(list=ls())
if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 2/")
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 2")
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 2/")
}
library(data.table)
#load data
data <- data.table(read.table("earinfect.txt",header=TRUE))

# Look at data
head(data)
summary(data)
data[,sum(persons), by = swimmer]
data[, sum(persons), by =location]
data[, sum(persons), by=age]
data[, sum(persons), by = sex]

data[persons == max(persons)]
data[,freq := infections/persons]
data[freq==max(freq)]
