rm(list=ls())
setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 2")
library(data.table)
library(ggplot2)
library(GGally)
#load data
data <- data.table(read.csv("clothing.csv",header=TRUE))

# Look at data
head(data)
tail(data)
summary(data)

# Remove index column
data$X <- NULL


ggpairs(data)
par(mfrow=c(2,3))
plot(clo~tOut,data=data,main="tOut")
plot(clo~tInOp,data=data,main="tInOp")
boxplot(clo~sex,data=data,main="sex")
boxplot(clo~subjId,data=data,main="subjId")
plot(clo~time,data=data,main="time")
boxplot(clo~day,data=data,main="day")

# Based on the plot:
#   - There seems to be a correlation between temperature and clothing
#   - Also seems to be effect of sex
#   - No effect from day or time or subjId

# Question: Does it make sense to use subjid???

# Let's look at y before choosing distribution
data$clo
# Values between 0 and 1 - kind of like percentage?
#     - Can't use binomial since it is counting number of successes
#     - According to the internet, use beta distribution
# Correction: according to ?binomial, the y-vector can be given as a numerical vector
# with values between 0 and 1 (proportion of successful cases). Though, then I have to give
# the total number of cases in weights (????)

data$weight <- 10
fit0 <- glm(clo ~ tOut + tInOp + sex + time + day, data = data, 
            family = binomial(link = "logit"), weights = data$weight)




