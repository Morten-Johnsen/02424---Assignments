rm(list=ls())
library(GGally)
library(qqplotr)
library(corrplot)
library(car)
library(reshape)
library(tidyverse)
library(magrittr)
library(dplyr)
library(betareg)
library(statmod)
library(jtools)
library(numDeriv)
library(latex2exp)
library(broom.mixed)
library(gridExtra)
library(grid)
library(lattice)
library(ggplot2)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 2/")
  figpath <- "/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 2/figs/"
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 2")
  figpath <- "~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 2/figs/"
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 2/")
  figpath <- "C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/Assignment 1/02424---Assignments/Assignment 2/"
}

earinfect <- read.table("earinfect.txt", header = T)
clothing <- read.csv(file = "clothing.csv", header = T)

#### Clothing ####
#1)
c.data <- dplyr::select(clothing, -subjId, -day, -time, -X)
head(c.data)

melt(c.data) %>%
  ggplot()+
  geom_histogram(aes(x = value, fill = sex), position = "identity", alpha = 0.4)+
  scale_fill_manual(values = c("blue", "orange"))+
  facet_wrap(~variable, scales = "free")+
  theme_bw()
ggsave(paste0(figpath,"clothing_data_histograms.png"), width = 30, height = 10, units = "cm")
#### 1) Distributions
beta.dist <- function(theta){
  return(-sum(dbeta(c.data$clo, shape1 = theta[1], shape2 = theta[2], log = T)))
}

norm.dist <- function(theta){
  return(-sum(dnorm(c.data$clo, mean = theta[1], sd = theta[2], log = T)))
}

gamma.dist <- function(theta){
  return(-sum(dgamma(c.data$clo, shape = theta[1], rate = theta[2], log = T)))
}

invgaus.dist <- function(theta){
  return(-sum(dinvgauss(c.data$clo, mean = theta[1], dispersion = theta[2], log =T)))
}

lnorm.dist <- function(theta){
  return(-sum(dlnorm(x=c.data$clo, meanlog = theta[1], sdlog = theta[2], log = T)))
}

beta.hat <- nlminb(start = c(1,1), objective = beta.dist)
norm.hat <- nlminb(start = c(1,1), objective = norm.dist)
gamma.hat <- nlminb(start = c(1,1), objective = gamma.dist)
invgaus.hat <- nlminb(start = c(1,1), objective = invgaus.dist)
lnorm.hat <- nlminb(start = c(1,1), objective = lnorm.dist)

ggplot(c.data)+
  geom_histogram(aes(x = clo, y = after_stat(density)), colour = "white", position = "identity", alpha = 0.4)+
  scale_fill_manual(values = c("blue", "orange"))+
  stat_function(aes(colour = "1: beta"), fun = dbeta, args = list(shape1 = beta.hat$par[1],shape2 = beta.hat$par[2]))+
  stat_function(aes(colour = "2: normal"), fun = dnorm, args = list(mean = norm.hat$par[1], sd = norm.hat$par[2]))+
  stat_function(aes(colour = "3: gamma"), fun = dgamma, args = list(shape = gamma.hat$par[1],rate = gamma.hat$par[2]))+
  stat_function(aes(colour = "4: gamma"), fun = dinvgauss, args = list(mean = invgaus.hat$par[1],dispersion = invgaus.hat$par[2]))+
  stat_function(aes(colour = "5: log-normal"), fun = dlnorm, args = list(meanlog = lnorm.hat$par[1],sdlog = lnorm.hat$par[2]))+
  theme_bw()+
  labs(y = "", colour = "Distribution")+
  scale_colour_manual(values = c("blue", "orange", "black", "red", "purple", "grey"), labels = c(paste0("Beta [AIC: ",round(2*beta.hat$objective + 2*length(beta.hat$par),3),"]"),
                                                                        paste0("Normal [AIC: ",round(2*norm.hat$objective + 2*length(norm.hat$par),3),"]"), 
                                                                        paste0("Gamma [AIC: ",round(2*gamma.hat$objective + 2*length(gamma.hat$par),3),"]"),
                                                                        paste0("Inv Gauss [AIC: ",round(2*invgaus.hat$objective + 2*length(invgaus.hat$par),3),"]"),
                                                                        paste0("Log-normal [AIC: ",round(2*lnorm.hat$objective + 2*length(lnorm.hat$par),3),"]")))+
  ggtitle("Clothing insulation level")
ggsave(paste0(figpath,"distribution_choice_1.png"), width = 20, height = 10, units = "cm")
#gamma is best (HUSK clo må gerne være > 1, det er ikke en fraktion).
fit.gamma <- glm(clo ~ tOut + tInOp + sex, data = c.data, family = Gamma(link = "cloglog"))
add1(object = fit.gamma, scope = ~.+tOut:sex+tInOp:sex+tOut:tInOp, test = "Chisq")
fit.gamma <- update(fit.gamma, .~.+tOut:sex)
Anova(fit.gamma, type = "III", test = "LR")
add1(object = fit.gamma, scope = ~.+tInOp:sex+tOut:tInOp, test = "Chisq")
drop1(object = fit.gamma, test = "Chisq")
fit.gamma <- update(fit.gamma, .~.-tInOp)
anova(fit.gamma, test = "Chisq")
Anova(fit.gamma, type = "III", test = "LR")
#goodness of fit:
pchisq(deviance(fit.gamma), df = dim(c.data)[1] - length(coefficients(fit.gamma)), lower.tail = F)
#er passende
summary(fit.gamma)
1 - pchisq(48.344,df=799) ## pass the goodness of fit.

confint(fit.gamma)

#Manuelt
glm.gamma.w <- function(theta){
  y <- c.data$clo
  k <- 1/theta[5]
  
  eta <- theta[1] + theta[2] * c.data$tOut + theta[3] * as.numeric(c.data$sex == "male") + theta[4] * as.numeric(c.data$sex == "male") * c.data$tOut
  mu <- 1-exp(-exp(eta))
  
  #using deviance
  #d <- 2*(y/mu - log(y/mu) - 1)
  #return(1/2 * sum(w*d))
  
  #using pdf and log-likelihood
  nll <- -sum(dgamma(y, shape = k, scale = mu/k, log = T))
  return(nll)
}

manual.fit <- nlminb(start = c(0,0,0,0,1), objective = glm.gamma.w)
manual.sd <- sqrt(diag(solve(hessian(func = glm.gamma.w, x = manual.fit$par))))
manual.fit$par
manual.sd

summary(fit.gamma)
#Det er mere eller mindre det samme, det er vist kun dispersion parameter der er lidt anderleedes
(AIC_manuel <- manual.fit$objective*2 + 5*2)
#manuel AIC er 0.03 bedre end glm.... det er det samme.

#### 2) residual analysis
par(mfrow=c(2,2))
postscript(file.path(getwd(), "fit.gamma.eps"), horizontal = FALSE, onefile = FALSE, paper = "special",height = 10, width = 10)
#png(filename = "/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 2/residual_analysis_1.png", width = 20, height = 10, units = "cm", res = 1000)
plot(fit.gamma, pch = 16)
dev.off()
c.data$residuals <- fit.gamma$residuals
c.data$leverage <- hatvalues(fit.gamma)
dev.off()

#gender-specific residual analysis
c.data$pred <- predict(fit.gamma)
c.data$pearson <- residuals(fit.gamma, type = "pearson")

sigma_sq <- fit.gamma$deviance / (dim(c.data)[1] - length(coefficients(fit.gamma)))
c.data$stdpearson <- c.data$pearson/sqrt(sigma_sq*(1-c.data$leverage))

first <- ggplot(c.data)+
  geom_point(aes(x = pred, y = residuals))+
  geom_hline(aes(yintercept = 0), colour = "blue", linetype = "dashed")+
  geom_smooth(aes(x = pred, y = residuals), colour = "blue", se = F)+
  theme_bw()+
  labs(x = "Predicted", y = "Residuals")+
  ggtitle("Residuals vs Fitted")

second <- ggplot(c.data)+
  geom_point(aes(x = pred, y = sqrt(stdpearson)))+
  geom_smooth(aes(x=pred,y = sqrt(stdpearson)), colour = "blue", se = F)+
  theme_bw()+
  labs(x = "Predicted", y = TeX("$\\sqrt{Std. Pearson Residuals}$"))+
  ggtitle("Scale-Location")

third <- ggplot(c.data, aes(sample = stdpearson))+
  stat_qq_band(fill = "blue", alpha = 0.2)+
  stat_qq_line(colour = "blue")+
  stat_qq_point()+
  theme_bw()+
  labs(x = "Theoretical quantiles", y = "Std. Pearson Residuals")+
  ggtitle("Normal QQ")

fourth <- ggplot(c.data, aes(x = leverage, y = stdpearson))+
  geom_point()+
  geom_hline(aes(yintercept = 0), colour = "blue", linetype = "dashed")+
  geom_smooth(colour = "blue", se = F)+
  theme_bw()+
  labs(x = "Leverage", y = "Std. Pearson Residuals")+
  ggtitle("Residuals vs Leverage")

grid.arrange(first, third, second, fourth, nrow = 2)

#Check whether residuals are i.i.d.
ggplot(c.data)+
  geom_boxplot(aes(x = sex, y = pearson, fill = sex))+
  theme_bw()+
  ggtitle("Residual variation difference between genders")
ggsave(filename = paste0(figpath,"gender_variance.png"), width = 20, height = 10, units = "cm")
#=> Residuals are not identically distributed. This will be addressed later.

par(mfrow=c(1,2))
p1 <- acf(c.data$pearson, main = "ACF pearson residuals")
p2 <- pacf(c.data$pearson, main = "PACF pearson residuals")
#=> residuals are not independent.

#### 3) Model interpretation
plot_summs(fit.gamma)+
  ggtitle("Parameter estimates [95% CI]")
ggsave(filename = paste0(figpath,"forest1.png"), width = 20, height = 10, units = "cm")

summary(fit.gamma)

#### 4) Fitting the model using subjId instead of sex
c.data2 <- dplyr::select(clothing, -X)
c.data2 %<>% 
  mutate(subjId = factor(subjId))
fit.gamma2 <- glm(clo ~ tOut + tInOp + subjId, data = c.data2, family = Gamma(link = "cloglog"))
anova(fit.gamma2, test = "Chisq")
Anova(fit.gamma2, type = "III")
add1(object = fit.gamma2, scope = ~.+I(tOut^2)+I(tInOp^2)+tOut:tInOp+tOut:subjId+tInOp:subjId, test = "Chisq")
fit.gamma2 <- update(fit.gamma2, .~.+tOut:subjId)
drop1(fit.gamma2, test = "Chisq")
anova(fit.gamma2, test = "Chisq")
# now all terms are significant and we see if additional terms should be added
add1(object = fit.gamma2, scope = ~.+I(tOut^2)+I(tInOp^2)+tOut:tInOp+tInOp:subjId, test = "Chisq")
fit.gamma2 <- update(fit.gamma2, .~.+tInOp:subjId)
Anova(fit.gamma2, type = "III")
drop1(fit.gamma2, test = "Chisq")
#type III anova show that tInOp is now no longer significant. We keep it, as it becomes significant again later
#all terms are now significant, see if we can add more information
add1(object = fit.gamma2, scope = ~.+I(tOut^2)+I(tOut^2):subjId+I(tInOp^2)+I(tInOp^2):subjId+tOut:tInOp+tOut:tInOp:subjId, test = "Chisq")
fit.gamma2 <- update(fit.gamma2, .~.+I(tInOp^2))
Anova(fit.gamma2, type = "III")
add1(object = fit.gamma2, scope = ~.+I(tOut^2)+I(tInOp^2)+tOut:tInOp+tOut:tInOp:subjId, test = "Chisq")
drop1(object = fit.gamma2, test = "Chisq")
#all terms are significant, no further terms can be added or removed.
Anova(fit.gamma2, type = "III")
fit.gamma2 <- update(fit.gamma2, .~.-tOut)
add1(object = fit.gamma2, scope = ~.+I(tOut^2)+I(tInOp^2)+tOut:tInOp+tOut:tInOp:subjId, test = "Chisq")
drop1(object = fit.gamma2, test = "Chisq")
summary(fit.gamma2)

#subject id is highly significant -> normally this would be an indicator to use a mixed model instead.
#the residual deviance for this model is a lot lower than for the sex-based model above.

1 - pchisq(12.715,df=661) ## pass the goodness of fit.

#### 5) Residual analysis including within day autocorrelation
c.data2$pearson <- residuals(fit.gamma2, type = "pearson")
ggplot(c.data2)+
  geom_boxplot(aes(x = subjId, y = pearson, fill = sex))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y = "Pearson residual", x = "Subject ID")
ggsave(paste0(figpath,"residual_subjid2.png"), width = 20, height = 10, units = "cm")

#still not constant residual variance as variations can be seen based on the subjId
c.data2$residuals <- fit.gamma2$residuals
c.data2$leverage <- hatvalues(fit.gamma2)
#gender-specific residual analysis
c.data2$pred <- predict(fit.gamma2)

sigma_sq <- fit.gamma2$deviance / (dim(c.data2)[1] - length(coefficients(fit.gamma2)))
c.data2$stdpearson <- c.data2$pearson/sqrt(sigma_sq*(1-c.data2$leverage))

first <- ggplot(c.data2)+
  geom_point(aes(x = pred, y = residuals))+
  geom_hline(aes(yintercept = 0), colour = "blue", linetype = "dashed")+
  geom_smooth(aes(x = pred, y = residuals), colour = "blue", se = F)+
  theme_bw()+
  labs(x = "Predicted", y = "Residuals")+
  ggtitle("Residuals vs Fitted")

second <- ggplot(c.data2)+
  geom_point(aes(x = pred, y = sqrt(stdpearson)))+
  geom_smooth(aes(x=pred,y = sqrt(stdpearson)), colour = "blue", se = F)+
  theme_bw()+
  labs(x = "Predicted", y = TeX("$\\sqrt{Std. Pearson Residuals}$"))+
  ggtitle("Scale-Location")

third <- ggplot(c.data2, aes(sample = stdpearson))+
  stat_qq_band(fill = "blue", alpha = 0.2)+
  stat_qq_line(colour = "blue")+
  stat_qq_point()+
  theme_bw()+
  labs(x = "Theoretical quantiles", y = "Std. Pearson Residuals")+
  ggtitle("Normal QQ")

fourth <- ggplot(c.data2, aes(x = leverage, y = stdpearson))+
  geom_point()+
  geom_hline(aes(yintercept = 0), colour = "blue", linetype = "dashed")+
  geom_smooth(colour = "blue", se = F)+
  theme_bw()+
  labs(x = "Leverage", y = "Std. Pearson Residuals")+
  ggtitle("Residuals vs Leverage")

png(filename = paste0(figpath,"residual_analysis_2.png"), width = 20, height = 10, units = "cm", res = 1000)
grid.arrange(first, third, second, fourth, nrow = 2)
dev.off()

#again the residuals are fairly well-behaved as compared to the predicted values
#showing a few outliers, but all residuals are within cooks distance.
#the qqplot show poor normality. Even when accounting for 95% CIs.
qqPlot(residuals(fit.gamma2))

#looking at the within day autocorrelation to examine independency.
par(mfrow=c(1,2))
acf(residuals(fit.gamma2, type = "pearson"), main = "ACF")
pacf(residuals(fit.gamma2, type = "pearson"), main = "PACF")
#still seeing overall autocorrelation
acf <- c()
lag <- c()
subject <- c()

for (i in unique(c.data2$subjId)){ # Loop over subjects
  # Take out data for that subject
  c.data2 %>% 
    filter(subjId == i) -> tmp
  for (j in unique(tmp$day)){ # Loop over different days that the subject has been at work
    # Look at that specific day, calculate acf
    tmp %>%
      filter(day == j) %>%
      arrange(time) %>%
      select(pearson) %>%
      acf(plot = F) -> acf.tmp
    # Save acf value for that subject on that day
    acf <- c(acf, acf.tmp$acf)
    lag <- c(lag, 0:(length(acf.tmp$acf)-1))
    subject <- c(subject, rep(i, length(acf.tmp$acf)))
  }
}

# Calculate the average per lag (accross subjects)
acf.data <- data.frame("acf" = acf, "lag" = lag, "subject" = subject)
mean.acf.data <- acf.data %>%
  group_by(lag) %>%
  summarise(lag_avg = mean(acf))
mean.acf.data

ggplot(acf.data)+
  geom_col(aes(x = lag, y = acf, fill = subject), position = "dodge", show.legend = F)+
  theme_bw()+
  scale_fill_manual(values = rep("grey",length(unique(subject))))+
  geom_hline(aes(yintercept = qnorm(1-0.05/2)/sqrt(dim(c.data2)[1]), colour = "95% significance level"), linetype = "dashed", linewidth = 0.8)+
  geom_hline(aes(yintercept = -qnorm(1-0.05/2)/sqrt(dim(c.data2)[1])), linetype = "dashed", colour = "royalblue1", linewidth = 0.8)+
  geom_segment(data = mean.acf.data, aes(x = lag-0.5, xend = lag+0.5, y = lag_avg, yend = lag_avg, colour = "ACF average pr lag"))+
  scale_x_continuous(n.breaks = 6)+
  labs(colour = "")+
  scale_colour_manual(values = c("royalblue1", "black"))+
  theme(legend.position = "top")+
  guides(colour = guide_legend(override.aes=list(linetype = c("dashed", "solid"), linewidth = c(0.5, 0.5))))+
  ggtitle("Within day autocorrelation for each subject")
ggsave(paste0(figpath,"within_day_autocor.png"), width = 20, height = 10, units = "cm")

#6) Optimal weight/dispersion parameter
library(caret)
library(stringr)
dummy <- dummyVars(" ~ .", data = c.data2)
new <- data.frame(predict(dummy, newdata = c.data2))

subjects <- length(names(new)[str_detect(names(new), pattern = "subjId")])-1

#design matrix

# X <- as.matrix(cbind(1, 
#            new$tInOp, 
#            new[, names(new)[str_detect(names(new), pattern = "subjId")][-1]],
#            new$tInOp^2,
#            new[, names(new)[str_detect(names(new), pattern = "subjId")]]*new$tOut,
#            new[, names(new)[str_detect(names(new), pattern = "subjId")][-1]]*new$tInOp))

X <- as.matrix(cbind(1, new$tOut, new$sexmale, new$tOut*new$sexmale))

glm.gamma.w2 <- function(theta){
  y <- new$clo
  #k <- rep(1/theta[1], length(y))
  beta <- theta[-c(1,2)]
  k <- numeric(length(y))
  
  k.male <- 1/theta[1]
  k.female <- 1/theta[2]
  k[as.logical(new$sexmale)] <- k.male
  k[as.logical(new$sexfemale)] <- k.female #shape
  
  eta <- as.numeric(X %*% matrix(beta, ncol = 1))

  mu <- 1-exp(-exp(eta))
  
  return(-sum(dgamma(y, shape = k, scale = mu/k, log = T)))
}

#manual.fit2 <- nlminb(start = c(0.02,as.numeric(coefficients(fit.gamma2))-0.001), objective = glm.gamma.w2, control=list(eval.max = 500, iter.max = 1000))
manual.fit2 <- nlminb(start = c(1,1,as.numeric(coefficients(fit.gamma))), objective = glm.gamma.w2)
#dispersion parameter males = 1/k_male, dispersion parameter females = 1/k_female, model parameters....
manual.fit2$par
1/manual.fit2$par[1:2]

manual.fit2$objective
summary(fit.gamma)

#looking at the variance associated with each gender based on the weight:
theta.hat <- c(manual.fit2$par[3:6], 1/manual.fit2$par[1:2])
#for the gamma distribution we get: var = shape*scale^2 = k*theta^2
#with our parametrization: scale = mu/shape = mu/k
#Thus: var = shape * mu^2/shape^2 = mu^2/shape = mu^2/k
#(our weight w is actually an estimate of the shape/dispersion parameter k) and thus:
##male
var.male <- (theta.hat[1] + theta.hat[2] * c.data$tOut[c.data$sex == "male"] + theta.hat[3] * as.numeric(c.data$sex == "male")[c.data$sex == "male"] + theta.hat[4] * as.numeric(c.data$sex == "male")[c.data$sex == "male"] * c.data$tOut[c.data$sex == "male"])^2 / theta.hat[5]
mean(var.male)
lines(var.male, col = 3)
##female
var.female <- (theta.hat[1] + theta.hat[2] * c.data$tOut[c.data$sex == "female"] + theta.hat[3] * as.numeric(c.data$sex == "female")[c.data$sex == "female"] + theta.hat[4] * as.numeric(c.data$sex == "female")[c.data$sex == "female"] * c.data$tOut[c.data$sex == "female"])^2 / theta.hat[6]
mean(var.female)
plot(var.female, type = "l")
#Higher variance for the women as expected

#Notes on the variance function V(.) as compareed to the variance Var(.) on page 93.

#7) Profile likelihood
glm.gamma.w.pf <- function(d_male,d_female){
  y <- new$clo
  k <- numeric(length(y))
  k.male <- 1/d_male
  k.female <- 1/d_female
  k[as.logical(new$sexmale)] <- k.male
  k[as.logical(new$sexfemale)] <- k.female #shape
  
  tmp.func <- function(beta, shape = k){
    eta <- as.numeric(X %*% matrix(beta, ncol = 1))
    
    mu <- 1-exp(-exp(eta))
    return(-sum(dgamma(y, shape = shape, scale = mu/shape, log = T)))
  }
  fit.tmp <- nlminb(start = c(0,0,0,0), objective = tmp.func)
  return(fit.tmp$objective)
}


w1 <- seq(0.025, 0.04, length.out = 20)
w2 <- seq(0.075, 0.095, length.out = 20)
z <- outer(w1, w2, FUN = Vectorize(function(w1,w2) glm.gamma.w.pf(w1,w2)))

png(filename = "countour_pf.png", width = 20, height = 10, units = "cm", res = 500)
contour(w1, w2, z, xlab = "male", ylab = "female", nlevels = 50)
dev.off()

plot_ly(x = w1, y = w2, z = z, type = "contour")
library(plotly)

glm.gamma.w.pf2 <- function(d){
  d_male <- d[1]
  d_female <- d[2]
  y <- new$clo
  k <- numeric(length(y))
  k.male <- 1/d_male
  k.female <- 1/d_female
  k[as.logical(new$sexmale)] <- k.male
  k[as.logical(new$sexfemale)] <- k.female #shape
  
  tmp.func <- function(beta, shape = k){
    eta <- as.numeric(X %*% matrix(beta, ncol = 1))
    
    mu <- 1-exp(-exp(eta))
    return(-sum(dgamma(y, shape = shape, scale = mu/shape, log = T)))
  }
  fit.tmp <- nlminb(start = c(0,0,0,0), objective = tmp.func)
  return(fit.tmp$objective)
}

z <- apply(reshape::expand(data.frame("d_male" = w1, "d_female" = w2)), MARGIN = 1, glm.gamma.w.pf2)

#### Earinfections ####
#1)
library(car)
head(earinfect)
earinfect$swimmer <- factor(earinfect$swimmer)
earinfect$location <- factor(earinfect$location)
earinfect$age <- factor(earinfect$age, ordered = T)
earinfect$sex <- factor(earinfect$sex)

fit.pois <- glm(infections ~ offset(persons) + age * sex * location * swimmer, data = earinfect, family = poisson(link = 'log'))
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-age:sex:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-age:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-age:sex:location:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-age:location:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-sex:location:swimmer)
anova(fit.pois, test = "Chisq")
fit.pois <- update(fit.pois, .~.-sex:location)
Anova(fit.pois, type = "III", test = "LR")
anova(fit.pois, test = "Chisq")


summary(fit.pois)
#check whether the model is appropriate or if there is overdispersion.
pchisq(fit.pois$deviance, df = fit.pois$df.residual, lower.tail = F) #model is appropriate
par(mfrow = c(2,2))
plot(fit.pois)
#for poisson glm with overdispersion see p. 152

#looking at binomial model:
earinfect$resp <- cbind(earinfect$infections, earinfect$persons - earinfect$infections)
fit.logis <- glm(resp ~ age + sex + location + swimmer, data = earinfect, family = binomial(link = "logit"))
anova(fit.logis, test = "Chisq")
Anova(fit.logis, type = "III", test.statistic = "LR")
fit.logis <- update(fit.logis, .~.-sex)
anova(fit.logis, test = "Chisq")
Anova(fit.logis, type = "III", test.statistic = "LR")
fit.logis <- update(fit.logis, .~.-age)
anova(fit.logis, test = "Chisq")
Anova(fit.logis, type = "III", test.statistic = "LR")
fit.logis <- update(fit.logis, .~.-swimmer)

add1(fit.logis, scope=~.+I(location:swimmer)+I(sex:location)+I(sex:age)+I(sex:swimmer)+I(age:swimmer)+I(location:swimmer:sex), test = "Chisq")

summary(fit.logis)
pchisq(fit.logis$deviance, df = fit.logis$df.residual, lower.tail = F) #model is highly appropriate
par(mfrow=c(2,2))
plot(fit.logis)
earinfect$logis.pred <- predict(fit.logis, type = "response")
ggplot(earinfect, aes(x = as.numeric(location)))+
  geom_point(aes(y = infections/persons))+
  geom_line(aes(y = logis.pred))+
  theme_bw()+
  scale_x_continuous(breaks = c(1,2), labels = c("beach", "non beach"))+
  ggtitle("Logistic regression model")



# Poisson family + added 'offset(log(persons))' ----------------------------------
# See page 123 in book

fit.pois <- glm(infections ~ offset(log(persons)) + age * sex * location * swimmer, data = earinfect, family = poisson)
fit.pois.full <- fit.pois
1-pchisq(3.3307e-15, 0)

drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-age:sex:location:swimmer)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-age:sex:location)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-age:location:swimmer)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-age:sex:swimmer)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-sex:location:swimmer)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-sex:location)
Anova(fit.pois, type = "III", test = "LR")
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-age:location)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-age:sex)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-sex:swimmer)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-age:swimmer)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois <- update(fit.pois, .~.-location:swimmer)
anova(fit.pois, test = "Chisq")
drop1(fit.pois)
fit.pois.final <- fit.pois

summary(fit.pois)

anova(fit.pois.full, fit.pois.final, test = "Chisq")
# We accept the reduced model because of the high p-value.

# Manual check
beta <- coefficients(fit.pois.final)
beta[is.na(beta)] <- 0
eta <- as.numeric(fit.pois.final$coefficients %*% beta + log(earinfect$infections))
all.equal(eta, predict(fit.pois.final))
# Soo not quite cool/good...

par(mfrow = c(2, 2))
plot(fit.pois.final)

1 - pchisq(fit.pois.final$deviance, df = fit.pois.final$df.residual, lower.tail = F) #model is highly appropriate

par(mfrow=c(1,2))
acf(residuals(fit.pois.final, type = "pearson"))
pacf(residuals(fit.pois.final, type = "pearson"))

# Log-distribution of persons
earinfect$logPersons <- log(earinfect$persons)
earinfect$loginfections <- log(earinfect$infections)
ggpairs(earinfect)


#### 1) Distributions
pois.dist <- function(theta){
  return(-sum(dpois(earinfect$infections, lambda = theta, log = T)))
}

neg.binom <- function(theta){
  return(-sum(dnbinom(earinfect$infections, size = theta[1], mu = theta[2], log = T)))
}

pois.hat <- nlminb(start = c(1), objective = pois.dist)
nbinom.hat <- nlminb(start = c(1,1), objective = neg.binom)
pois.dens <- dpois(1:21, lambda = pois.hat$par)
nbinom.dens <- dnbinom(1:21, size = nbinom.hat$par[1], mu = nbinom.hat$par[2])

ggplot(earinfect)+
  geom_histogram(aes(x = infections, y = after_stat(density)), colour = "white", position = "identity", alpha = 0.4, bins = 9)+
  scale_fill_manual(values = c("blue", "orange"))+
  geom_step(data = data.frame(pois.dens), aes(x = 1:21, y = pois.dens, colour = "1: poisson"))+
  geom_step(data = data.frame(nbinom.dens), aes(x = 1:21, y = nbinom.dens, colour = "2: negative binomial"))+
  theme_bw()+
  theme(legend.position = "top")+
  labs(y = "", colour = "Distribution", x = "Ear infections")+
  scale_colour_manual(values = c("blue", "orange", "black", "red", "purple", "grey"), 
                      labels = c(paste0("Poisson [AIC: ",round(2*pois.hat$objective + 2*length(pois.hat$par),3),"]"),
                                 paste0("Negative Binomial [AIC: ",round(2*nbinom.hat$objective + 2*length(nbinom.hat$par),3),"]")))+
  ggtitle("Ear infections")
ggsave("distribution_choice_3.png", width = 20, height = 10, units = "cm")



# Plot different disitrbutions --------------------------------------------

library(fitdistrplus)
par(mfrow = c(1, 1))
fit.lnorm <- fitdist(earinfect$infections, "lnorm")
fit.negbin <- fitdist(earinfect$infections, "nbinom")
fit.pois <- fitdist(earinfect$infections, "pois")
fit.gamma <- fitdist(earinfect$infections, "gamma")
fit.exp <- fitdist(earinfect$infections, "exp")
plot.legend <- c("Log-normal","Negative binomial", "Poisson", "Gamma", "Exponential")
denscomp(list(fit.lnorm, fit.negbin, fit.pois, fit.gamma, fit.exp), 
         xlab = "Earinfections",fitlwd = 2, legendtext = plot.legend, 
         xlim = c(0, max(earinfect$infections)), ylim = c(0, 0.2))
qqcomp(list(fit.lnorm, fit.negbin, fit.pois, fit.gamma), legendtext = plot.legend, xlim = c(0, max(earinfect$infections)))

cat("AIC (log-normal):", fit.lnorm$aic, 
    "\nAIC (negative binomial):", fit.negbin$aic,
    "\nAIC (poisson):", fit.pois$aic, 
    "\nAIC (gamma):", fit.gamma$aic, 
    "\nAIC (exponential):", fit.exp$aic)


# plots -------------------------------------------------------------------

ggplot(earinfect, aes(x = infections, fill = factor(location))) +
  geom_bar(stat = "count") +
  scale_fill_discrete(name = "Location",
                      labels = c("Beach", "Non Beach"))

ggplot(earinfect, aes(x = infections, fill = factor(swimmer))) +
  geom_bar(stat = "count") +
  scale_fill_discrete(name = "Swimmer",
                      labels = c("Frequent", "Occational"))

ggplot(earinfect, aes(x = infections, fill = factor(age))) +
  geom_bar(stat = "count") +
  scale_fill_discrete(name = "Age groups",
                      labels = c( "15-19", "20-24", "25-29"))

ggplot(earinfect, aes(x = infections, fill = factor(sex))) +
  geom_bar(stat = "count") +
  scale_fill_discrete(name = "Sex",
                      labels = c( "Female", "Male"))



# offset log or not log ---------------------------------------------------


plot(earinfect$infections ~ earinfect$persons)
abline(lm(infections ~ persons,data=earinfect),col='red')

par(mfrow = c(1, 2))
plot(earinfect$infections ~ log(earinfect$persons))
abline(lm(infections ~ log(persons),data=earinfect),col='red')

summary(lm(infections ~ log(persons),data=earinfect))$r.squared
summary(lm(infections ~ persons,data=earinfect))$r.squared




# negative binomial -------------------------------------------------------

fitnb <- glm.nb(infections ~ offset(persons) + age * sex * location * swimmer, 
                data = earinfect)
summary(fitnb)
drop1(fitnb)
fitnb <- update(fitnb, .~.-age:sex:location:swimmer)
summary(fitnb)
drop1(fitnb)
fitnb <- update(fitnb, .~.-age:sex:swimmer)
summary(fitnb)
drop1(fitnb)
fitnb <- update(fitnb, .~.-age:location:swimmer)
summary(fitnb)
drop1(fitnb)
fitnb <- update(fitnb, .~.-age:swimmer)
summary(fitnb)
drop1(fitnb)
fitnb <- update(fitnb, .~.-location:swimmer)
summary(fitnb)
drop1(fitnb)
fitnb <- update(fitnb, .~.-sex:location:swimmer)
summary(fitnb)
drop1(fitnb)
fitnb <- update(fitnb, .~.-sex:swimmer)
summary(fitnb)
drop1(fitnb)
fitnb <- update(fitnb, .~.-swimmer)
summary(fitnb)

par(mfrow = c(2, 2))
plot(fitnb)

# Hmmmm


dummy <- dummyVars(" ~ .", data = earinfect)
earinf <- data.frame(predict(dummy, newdata = earinfect))

X <- cbind(1, earinf$swimmerOccas, earinf$locationNonBeach, earinf$sexMale, earinf$age20.24, earinf$age25.29)

nbinom.glm <- function(theta){
  trials <- earinf$persons
  successes <- earinf$infections
  p <- theta[1]
  beta <- theta[-1]
  eta <- as.numeric(X %*% matrix(beta, ncol = 1))
  
  mu <- exp(eta)
  
  return(-sum(dnbinom(x = successes, size = mu*(1-p)/(p), prob = 1-p, log = T)))
}

fit <- nlminb(start = c(0.5,rep(0.1,6)), objective = nbinom.glm)

library(MASS)
fit2 <- glm.nb(infections ~ swimmer + location + sex + age, 
               link = log, 
               offset = persons,
               data = earinfect)
fit$par
fit2$coefficients
summary(fit2)
