rm(list=ls())
setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 2")
library(data.table)
#load data
data <- data.table(read.csv("clothing.csv",header=TRUE))

# Look at data
head(data)
summary(data)


