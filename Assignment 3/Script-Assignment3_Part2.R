######################## Part 2 ########################
rm(list=ls())

library(reshape)
library(tidyverse)
library(magrittr)
library(dplyr)
library(ggplot2)
library(nlme)
library(lme4)

if (Sys.getenv('USER') == "mortenjohnsen"){
  setwd("/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 3/")
  figpath <- "/Users/mortenjohnsen/OneDrive - Danmarks Tekniske Universitet/DTU/10. Semester/02424 - Advanced Dataanalysis and Statistical Modellling/02424---Assignments/Assignment 3/figs/"
}else if (Sys.getenv('USER') == "freja"){
  setwd("~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 3")
  figpath <- "~/Documents/Uni/TiendeSemester/Adv. data analysis and stat. modelling/02424---Assignments/Assignment 3/figs/"
}else{
  setwd("C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/02424---Assignments/Assignment 3/")
  figpath <- "C:/Users/catdu/OneDrive/DTU/10. semester/Advanced Dataanalysis and Statistical Modelling/02424---Assignments/Assignment 3/figs"
}

clothing <- read.csv(file = "clothingFullAss03.csv", header = T)
head(clothing)

unique(clothing$day)
unique(clothing$time)
c.data <- dplyr::select(clothing, -time2, -X)
head(c.data)

c.data$subjId <- factor(c.data$subjId)
c.data$day <- factor(c.data$day)
c.data$subDay <- factor(c.data$subDay)
c.data$sex <- factor(c.data$sex)

#### Part 2.1 ####
# Fitting model using lmer
fit0 <- lmer(clo~sex+(1|subjId),data=c.data,REML=FALSE)
summary(fit0)

# Create design matrices
X <- cbind(1, as.numeric(c.data$sex == "male"))
Z <- dummy(c.data$subjId, levelsToKeep = unique(c.data$subjId))
y <- c.data$clo
dim(Z) # Z is 803 by 47. 
# Row corresponds to observation for a certain subject on a certain time.

## Joint Likelihood
library(mvtnorm)
library(numDeriv)

# Function that outputs negative log-likelihood of multivariate normal distribution
opt.fun <- function(theta){
  # Create variance-covariance structures
  Psi <- diag(rep(exp(theta[1]),dim(Z)[2])) # random effects
  Sigma <- diag(rep(exp(theta[2]),length(y))) # fixed effects
  
  # fixed effect params
  beta <- matrix(theta[3:4], ncol = 1)
  
  # Dispersion matrix of clo (eq 13 in report, eq. 5.48 in book)
  V <- Sigma + Z%*%Psi%*%t(Z)
  
  #log-likelihood of multivariate normal distribution:
  obj <-  mvtnorm::dmvnorm(y, mean = X%*%beta, sigma = V, log = T)
  
  return(-obj)
}

# Optimize sigma_u^2, sigma^2, beta
# Right now we start with the values from the lmer-function to be close to the 
# true values
# The inputs to sigma_u^2 and sigma^2 are in the log-domain (since they are negative)
par_est <- nlminb(start = c(log(0.11665),log(0.09883), 0.59176, -0.08322),
                  objective = opt.fun, control = list(trace = 1))
# Get values of sigma_u and sigma in original domain
sqrt(exp(par_est$par[1:2]))
sds <- sqrt(diag(solve(hessian(func = opt.fun, x = par_est$par))))

# Compare to true values
manual_fit0 <- data.frame("Parameter" = c("Psi (Sigma.u)", "Sigma", "beta0", "beta1"),
                          "Interpretation" = c("SubjId", "Model Residual", "Model intercept", "Model slope (sex)"),
                          "Estimate" = c(exp(par_est$par[1:2]), par_est$par[3:4]))
cat("Log-likelihood = ", -par_est$objective)
manual_fit0
summary(fit0)

Psi <- diag(rep(exp(par_est$par[1]),dim(Z)[2]))
Sigma <- diag(rep(exp(par_est$par[2]),length(y)))
beta <- matrix(par_est$par[3:4], ncol = 1)


# Estimate random effects based on the found Psi, Sigma and beta
y_adj <- y - X%*%beta
# Eq. 5.60 in book
u <- solve(t(Z)%*%solve(Sigma)%*%Z + solve(Psi)) %*% (t(Z)%*%solve(Sigma)%*%y_adj)
# Compare to fit
cbind(ranef(fit0)$subjId,u)
#All values are the same :)))

#### 2.2 ####
# Fit nested model using lmer
fit1 <- lmer(clo~sex+(1|subjId)+(1|subjId:day),data=c.data,REML=FALSE)
ranef(fit1)

# Design matrices
X <- cbind(1, as.numeric(c.data$sex == "male"))
# Z is the same as before
Z <- dummy(c.data$subjId, levelsToKeep = unique(c.data$subjId))
# Introduce W which is the design matrix for the nested subjdID:day observations
c.data$subjDay <- factor(paste0(c.data$subjId, "/",c.data$day))
W <- dummy(c.data$subjDay, levelsToKeep = unique(c.data$subjDay))
y <- c.data$clo
dim(W) # W is 803 by 136 
# 136 = subjId*day

# theta = [log(sigma_v^2), log(sigma_u^2), log(sigma^2), beta_intercept, beta_slope ]
opt.fun1 <- function(theta){
  # Variance-Covariance matrices
  Psi <- diag(rep(exp(theta[2]),dim(Z)[2]))
  Sigma <- diag(rep(exp(theta[3]),length(y)))
  Phi <- diag(rep(exp(theta[1]),dim(W)[2]))
  
  beta <- matrix(theta[4:5], ncol = 1)
  
  # Same as before, but with + W%*%Phi%*%t(W) for the random effs. on subjID/Day
  V <- Sigma + Z%*%Psi%*%t(Z) + W%*%Phi%*%t(W)
  
  obj <-  mvtnorm::dmvnorm(y, mean = X%*%beta, sigma = V, log = T)
  
  return(-obj)
}

# Optimize params
par_est1 <- nlminb(start = c(0, 0, 0, 1, 1),
                   objective = opt.fun1, control = list(trace = 1))
# Compare to fit1
sds1 <- sqrt(diag(solve(hessian(func = opt.fun1, x = par_est1$par))))
summary(fit1)
manual_fit1 <- data.frame("Parameter" = c("Phi (Sigma.v)", "Psi (Sigma.u)", "Sigma", "beta0", "beta1"),
                          "Interpretation" = c("SubjId/day", "SubjId", "Model Residual", "Model Intercept", "Model Slope (sex)"),
                          "Estimate" = c(exp(par_est1$par[1:3]), par_est1$par[4:5]))
cat("Log-likelihood = ", -par_est1$objective)
manual_fit1
sds1

# Hierarchical estimation: We first estimated Beta, Psi and Phi.
# now we will estimate u and the nested value v.
Phi <- diag(rep(exp(par_est1$par[1]),dim(W)[2]))
Psi <- diag(rep(exp(par_est1$par[2]),dim(Z)[2]))
Sigma <- diag(rep(exp(par_est1$par[3]),length(y)))
beta <- matrix(par_est1$par[4:5], ncol = 1)

# Use equation on the bottom of page 184 in book.
A <- rbind(
  cbind(t(Z)%*%solve(Sigma)%*%Z+solve(Psi), t(Z)%*%solve(Sigma)%*%W),
  cbind(t(W)%*%solve(Sigma)%*%Z, t(W)%*%solve(Sigma)%*%W+solve(Phi))
)
B <- rbind(t(Z)%*%solve(Sigma)%*%(y-X%*%beta),
           t(W)%*%solve(Sigma)%*%(y-X%*%beta))
# Solve: A[1,1]*u + A[1,2]*v = B[1] for u
# Solve: A[2,2]*v + A[2,1]*u = B[2] for v
uv <- solve(A, B)

u <- t(t(uv[1:dim(Z)[2]]))
v <- t(t(uv[-c(1:dim(Z)[2]),]))

# Compare with fit: all values are correct :)))
cbind(ranef(fit1)$subjId, u)
cbind(ranef(fit1)$`subjId:day`, v)

data.frame(uv) %>%
  rownames_to_column(var = "SubDay") %>%
  mutate(uv = signif(uv, digits = 4),
         SubjectID = case_when(str_detect(SubDay, "/") ~ str_remove(SubDay, "/."), TRUE ~ SubDay),
         day = case_when(str_detect(SubDay, "/") & str_length(SubjectID) == 1 ~ paste0("Day ", str_remove(SubDay, "./"), " [v]"),
                         str_detect(SubDay, "/") & str_length(SubjectID) == 2 ~ paste0("Day ", str_remove(SubDay, "../"), " [v]"),
                         TRUE ~ "Subject [u]")) %>%
  remove_rownames() %>%
  select(uv, SubjectID, day) %>%
  pivot_wider(id_cols = SubjectID,
              names_from = day,
              values_from = uv) -> appendix_table_uv
library(xtable)
print(xtable(appendix_table_uv, display = c("d", "d", "g", "g", "g", "g", "g")
             , type = "latex"
             , digits = 4
             , caption = "Outliers excluded from the data set."
             , label = "tab:uv_table"),
      file = "uv_table.tex",
      caption.placement = "top", include.rownames = F, hline.after = c(-1,0,0,47))

data.frame("Residuals" = residuals(fit1, "pearson"), "Sex" = c.data$sex, "SubjID" = c.data$subjId) %>%
  ggplot()+
  geom_boxplot(aes(x = SubjID, y = Residuals, fill = Sex))+
  theme_bw()+
  labs(y = "Pearson residuals", x = "Subject ID")+
  ggtitle("Gender specific residual variance")
ggsave("compareVarianceSex.png", width = 20, height = 10, units = "cm")



#### 2.3 ####
#Implementation of sex specific variances.
#Reasonable implementations: Domain? Which parameters?
#Obviously estimation will be performed in the log-domain to avoid negative values
#Parameters to be identified: Scaling factor based on gender - Same scaling factor for all variances?
#Scaling factor alpha(sex): Assume there is a difference in variance across genders
#e.g. clothing insulation level vary less/more for women than for men.

#get the sex of each subj
c.data %>%
  select(subjId, sex) %>%
  distinct() %>%
  mutate(sex = case_when(sex == "male" ~ 1, TRUE ~ 0)) -> subjSex

c.data %>%
  select(subjId, day, sex) %>%
  distinct() %>%
  mutate(sex = case_when(sex == "male" ~ 1, TRUE ~ 0)) -> subjDaySex

# When is the observation a male
idx_Psi <- subjSex$sex
idx_Phi <- subjDaySex$sex
idx_Sigma <- as.numeric(c.data$sex == "male")

# theta = [alpha, log(sigma_v^2), log(sigma_u^2), log(sigma^2), beta_intercept, beta_slope ]
opt.fun2 <- function(theta){
  # alpha for males
  alpha <- theta[1] #constant effect of gender across all variances
  
  #Estimation in the exponential domain and estimate alpha as a scaling factor of the female variance
  # -> e(theta)*e(alpha) = e(theta + alpha)
  Psi <- diag(exp(theta[3] + idx_Psi*alpha))
  Sigma <- diag(exp(theta[4] + idx_Sigma*alpha))
  Phi <- diag(exp(theta[2] + idx_Phi*alpha))
  
  beta <- matrix(theta[5:6], ncol = 1)
  
  # Dispersion of clo
  V <- Sigma + Z%*%Psi%*%t(Z) + W%*%Phi%*%t(W)
  
  # Margial log-likelihood
  obj <-  mvtnorm::dmvnorm(y, mean = X%*%beta, sigma = V, log = T)
  
  return(-obj)
}

# Optimize params
par_est2 <- nlminb(start = c(0, 0, 0, 0, 1, 1),
                   objective = opt.fun2, control = list(trace = 1))
# Standard deviations of parameters
sds2 <- sqrt(diag(solve(hessian(opt.fun2, x = par_est2$par))))
# Look at values
manual_fit2 <- data.frame("Parameter" = c("alpha", "Phi (Sigma.v)", "Psi (Sigma.u)", "Sigma", "beta0", "beta1"),
                          "Interpretation" = c("Male gender variance scaling", "SubjId/day", "SubjId", "Model Residual", "Model Intercept", "Model Slope (sex)"),
                          "Estimate" = c(exp(par_est2$par[1:4]), par_est2$par[5:6]))
cat("Log-likelihood = ", -par_est2$objective)
manual_fit2
sds2

alpha <- par_est2$par[1] #constant effect of gender across all variances

#Estimation in the exponential domain and estimate alpha as a scaling factor of the female variance
# -> e(theta)*e(alpha) = e(theta + alpha)
beta <- matrix(par_est2$par[5:6], ncol = 1)
Psi <- diag(exp(par_est2$par[3] + idx_Psi*alpha))
Sigma <- diag(exp(par_est2$par[4] + idx_Sigma*alpha))
Phi <- diag(exp(par_est2$par[2] + idx_Phi*alpha))

# Find random effects
A2 <- rbind(
  cbind(t(Z)%*%solve(Sigma)%*%Z+solve(Psi), t(Z)%*%solve(Sigma)%*%W),
  cbind(t(W)%*%solve(Sigma)%*%Z, t(W)%*%solve(Sigma)%*%W+solve(Phi))
)
B2 <- rbind(t(Z)%*%solve(Sigma)%*%(y-X%*%beta),
           t(W)%*%solve(Sigma)%*%(y-X%*%beta))
uv2 <- solve(A2, B2)

data.frame(uv2) %>%
  rownames_to_column(var = "SubDay") %>%
  mutate(uv = signif(uv, digits = 4),
         SubjectID = case_when(str_detect(SubDay, "/") ~ str_remove(SubDay, "/."), TRUE ~ SubDay),
         day = case_when(str_detect(SubDay, "/") & str_length(SubjectID) == 1 ~ paste0("Day ", str_remove(SubDay, "./"), " [v]"),
                         str_detect(SubDay, "/") & str_length(SubjectID) == 2 ~ paste0("Day ", str_remove(SubDay, "../"), " [v]"),
                         TRUE ~ "Subject [u]")) %>%
  remove_rownames() %>%
  select(uv, SubjectID, day) %>%
  pivot_wider(id_cols = SubjectID,
              names_from = day,
              values_from = uv) -> appendix_table_uv2
library(xtable)
print(xtable(appendix_table_uv2, display = c("d", "d", "g", "g", "g", "g", "g")
             , type = "latex"
             , digits = 4
             , caption = "Estimates of the latent variables, u and v, when using a gender-specific variance."
             , label = "tab:uv_table2"),
      file = "uv_table2.tex",
      caption.placement = "top", include.rownames = F, hline.after = c(-1,0,0,47))

uv_opt0 <- uv2
#### 2.4 ####
#derivation in overleaf

#### 2.5 ####
c.data %>% 
  select(subjId, sex) %>%
  distinct() %>%
  mutate(sex = case_when(sex == "male" ~ 1, TRUE ~ 0)) -> subjSex

c.data %>% 
  select(subjId, day, sex) %>%
  distinct() %>%
  mutate(sex = case_when(sex == "male" ~ 1, TRUE ~ 0)) -> subjDaySex

# Number of observations per subDay
c.data %>%
  select(subjId, day) %>%
  distinct() %>%
  group_by(subjId) %>%
  summarise(subjDay_reps = sum(n())) %>%
  select(subjDay_reps) %>%
  t() %>% as.numeric() -> gamma_idxSubDay

# Number of observations per subjId
c.data %>%
  select(subjId) %>%
  group_by(subjId) %>%
  summarise(subj_reps = sum(n())) %>%
  select(subj_reps) %>%
  t() %>% as.numeric() -> gamma_idxSub

idx_Psi <- subjSex$sex
idx_Phi <- subjDaySex$sex
idx_Sigma <- as.numeric(c.data$sex == "male")

# Design matrices (same as 2.3)
X <- cbind(1, as.numeric(c.data$sex == "male"))
Z <- dummy(c.data$subjId, levelsToKeep = unique(c.data$subjId))
#Introduce W which is the design matrix for the nested subjdID:day observations
c.data$subjDay <- factor(paste0(c.data$subjId, "/",c.data$day))
W <- dummy(c.data$subjDay, levelsToKeep = unique(c.data$subjDay))
y <- c.data$clo


X <- cbind(1, as.numeric(c.data$sex == "male"))
Z <- dummy(c.data$subjId, levelsToKeep = unique(c.data$subjId))
#Introduce W which is the design matrix for the nested subjdID:day observations
c.data$subjDay <- factor(paste0(c.data$subjId, "/",c.data$day))
W <- dummy(c.data$subjDay, levelsToKeep = unique(c.data$subjDay))
y <- c.data$clo


joined.likelihood <- function(gamma, phi, Phi, Psi, Sigma, X, Z, W, beta, y){
  # Scale covariance-matriced with exp(-gamma) as in model formulation
  Phi <- Phi/exp(gamma)
  Psi <- Psi/exp(gamma)
  Sigma <- Sigma/exp(gamma)
  
  # Dispersion for clo
  V <- Sigma + Z%*%Psi%*%t(Z) + W%*%Phi%*%t(W)
  
  #constant effect of gender across all variances
  #Estimation in the exponential domain and estimate alpha as a scaling factor of the female variance
  # -> e(theta)*e(alpha) = e(theta + alpha)
  # Return negative joint likelihood: 
  #     -f_(clo,gamma) = -(f(clo|gamma)*f(gamma))
  # ->  -log(f_(clo,gamma)) = -(log(f(clo|gamma)) + log(f(gamma)))
  
  return(-(mvtnorm::dmvnorm(x = y, mean = X%*%beta, sigma = V, log = T)+dgamma(exp(gamma), shape = phi, rate = phi, log = T)))
}

# Phi <- diag(exp(theta[2] + idx_Phi*alpha + rep(-gamma, times = gamma_idxSubDay)))
# Psi <- diag(exp(theta[3] + idx_Psi*alpha - gamma))
# Sigma <- diag(exp(theta[4] + idx_Sigma*alpha + rep(-gamma, times = gamma_idxSub)))

opt.fun4 <- function(theta, X, Z, W, y){
  # alpha for males
  alpha <- theta[1] #constant effect of gender across all variances
  
  #Estimation in the exponential domain and estimate alpha as a scaling factor of the female variance
  # -> e(theta)*e(alpha) = e(theta + alpha)
  Phi <- diag(exp(theta[2] + idx_Phi*alpha))
  Psi <- diag(exp(theta[3] + idx_Psi*alpha))
  Sigma <- diag(exp(theta[4] + idx_Sigma*alpha))
  
  beta <- matrix(theta[5:6], ncol = 1)
  phi <- 1+exp(theta[7])
  nu <- phi * 2
  
  ##### OPTIMIZE GAMMA ####
  gAmmA <- c()
  #find gamma:
  # Loop over subject IDs to solve smaller problems, since they are independent.
  for (i in 0:(length(unique(c.data$subjId))-1)){
    idx <- which(c.data$subjId == i)
    # Days for this subject
    subday <- which(unique(c.data$subjDay) %in% as.character(unique(c.data$subjDay[idx])))
    days <- length(subday)
    
    # Estimate optimal gamma
    est <- nlminb(start = 0, objective = joined.likelihood,
                  #de her matricer er mere eller mindre de samme for alle iteration. Da de alle sammen bare er
                  #diagonal matricer og det eneste der ændrer sig er hvis der er !=18 observationer for en person
                  #eller hvis der er != 3 forskellige dage med observationer
                  
                  # Phi: covariances for subDay 
                  #   subset: Phi[subday,subday]
                  # X: Take rows and columns corresponding to subject
                  # Z: Take rows and column corresponding to subject
                  # W: Take rows and columns corresponding to subject
                  
                  phi = phi, Phi = Phi[subday,subday], Psi = matrix(Psi[i+1,i+1]), Sigma =Sigma[idx,idx],
                  X = X[idx,], Z = matrix(Z[idx,(i+1)], ncol = 1), W = W[idx, subday],
                  beta = beta, y = y[idx])
    gAmmA <- c(gAmmA, exp(est$par))
  }
  
  #########################
  
  # Dispersion of clo
  Phi <- Phi*diag(rep(1/gAmmA, times = gamma_idxSubDay))
  Psi <- Psi*diag(1/gAmmA)
  Sigma <- Sigma*diag(rep(1/gAmmA, times = gamma_idxSub))
  
  V <- Sigma + Z%*%Psi%*%t(Z) + W%*%Phi%*%t(W)
  
  # Marginal log-likelihood
  mu <- X%*%beta

  obj <-  mvtnorm::dmvt(x = as.numeric(y - mu), sigma = V*(nu - 2)/nu, df = nu, log = T)
  #obj <-  mvtnorm::dmvnorm(y - X%*%beta, sigma = V, log = T)
  cat("--------------------\n")
  cat(obj,"\n")
  cat(own.dmvt(x = y-mu, Sigma = V*(nu - 2)/nu, nu = nu),"\n")
  cat("--------------------\n")
  # 
  return(-obj)
}

opt.fun4(par_est4$par, X = X, Z = Z, W = W, y = y)
library(msos)
own.dmvt <- function(x, Sigma, nu){
  x <- matrix(x, nrow = length(x))
  p <- nrow(x)
  #return(gamma(1/2*(df+p)) * (1 + 1/df*t(x)%*%solve(Sigma)%*%(x))^(-(df+p)/2) / ((df*pi)^(p/2)*sqrt(det(Sigma))*gamma(df/2)))
  
  #return(lgamma((nu+p)/2) - lgamma(nu/2) - p/2*log(nu*pi) - 1/2*logdet(Sigma) - (nu+p)/2 * log(1+1/nu*t(x)%*%solve(Sigma)%*%x))
  return(lgamma((nu+p)/2) - (lgamma(nu/2) + p/2*log(nu*pi) + 1/2*logdet(Sigma)) - (nu+p)/2 * log(1+1/nu*(t(x)%*%solve(Sigma)%*%x)))
}
dec <- chol(Sigma)
lgamma((p + nu)/2) - (lgamma(nu/2) + sum(log(diag(dec))) + p/2 * log(pi * nu)) - 0.5 * (nu + p) * log1p(rss/nu)

mvtnorm::dmvt(x = as.numeric(y - mu), sigma = V*(nu - 2)/nu, df = nu, log = T)

# Optimize params
#Parameter rækkefølge: alpha, exp(sigma.v^2), exp(sigma.u^2), exp(sigma^2), beta0, beta1, 1+exp(phi)
par_est4 <- nlminb(start = rep(0, 7),
                   objective = opt.fun4, 
                   X = X, Z = Z, W = W, y = y,
                   control = list(trace = 1))
par_est4
# Standard deviations of parameters
sds4 <- sqrt(diag(solve(hessian(opt.fun4, x = par_est4$par, X = X, Z = Z, W = W, y = y))))

#### 2.6 ####
# Objective function to estimate gamma for a certain subject ID
# Inputs are subsets corresponding to the subject
joined.likelihood <- function(gamma, sigma.g, Phi, Psi, Sigma, X, Z, W, beta, y){
  # Scale covariance-matriced with exp(-gamma) as in model formulation
  Phi <- Phi*exp(-gamma)
  Psi <- Psi*exp(-gamma)
  Sigma <- Sigma*exp(-gamma)
  
  # Dispersion for clo
  V <- Sigma + Z%*%Psi%*%t(Z) + W%*%Phi%*%t(W)
  
  #constant effect of gender across all variances
  #Estimation in the exponential domain and estimate alpha as a scaling factor of the female variance
  # -> e(theta)*e(alpha) = e(theta + alpha)
  # Return negative joint likelihood: 
  #     -f_(clo,gamma) = -(f(clo|gamma)*f(gamma))
  # ->  -log(f_(clo,gamma)) = -(log(f(clo|gamma)) + log(f(gamma)))
  
  return(-(mvtnorm::dmvnorm(y, mean = X%*%beta, sigma = V, log = T)+dnorm(gamma, sd = sqrt(sigma.g), log = T)))
}

# This function is used after finding all gammas (one for each subject)
# It returns the negative joint likelihood
joined.likelihood.H <- function(gamma, sigma.g, Phi, Psi, Sigma, X, Z, W, beta, y){
  # Scale with exp(-gamma)
  Phi <- Phi*diag(exp(rep(-gamma, times = gamma_idxSubDay)))
  Psi <- Psi*diag(exp(-gamma))
  Sigma <- Sigma*diag(exp(rep(-gamma, times = gamma_idxSub)))
  
  # Dispersion matrix of y
  V <- Sigma + Z%*%Psi%*%t(Z) + W%*%Phi%*%t(W)
  
  # Return negative joint likelihood: 
  #     -f_(clo,gamma) = -(f(clo|gamma)*f(gamma))
  # ->  -log(f_(clo,gamma)) = -(log(f(clo|gamma)) + log(f(gamma)))
  return(-(mvtnorm::dmvnorm(y, mean = X%*%beta, sigma = V, log = T)+sum(dnorm(gamma, sd = sqrt(sigma.g), log = T))))
}

# theta = [alpha, log(sigma_v^2), log(sigma_u^2), log(sigma^2), log(sigma_g^2), 
#           beta_intercept, beta_slope]
opt.fun3 <- function(theta, X, Z, W, y){
  alpha <- theta[1]
  Phi <- diag(exp(theta[2] + idx_Phi*alpha))
  Psi <- diag(exp(theta[3] + idx_Psi*alpha))
  Sigma <- diag(exp(theta[4] + idx_Sigma*alpha))
  sigma.g <- exp(theta[5])
  beta <- matrix(theta[6:7], ncol = 1)
  gamma <- c()
  H <- c()
  
  # Do laplace approximation of marginal likelihood of clo
  
  #find gamma:
  # Loop over subject IDs to solve smaller problems, since they are independent.
  for (i in 0:(length(unique(c.data$subjId))-1)){
    idx <- which(c.data$subjId == i)
    # Days for this subject
    subday <- which(unique(c.data$subjDay) %in% as.character(unique(c.data$subjDay[idx])))
    days <- length(subday)
    
    # Estimate optimal gamma
    est <- nlminb(start = 0, objective = joined.likelihood,
                  #de her matricer er mere eller mindre de samme for alle iteration. Da de alle sammen bare er
                  #diagonal matricer og det eneste der ændrer sig er hvis der er !=18 observationer for en person
                  #eller hvis der er != 3 forskellige dage med observationer
                  
                  # Phi: covariances for subDay 
                  #   subset: Phi[subday,subday]
                  # X: Take rows and columns corresponding to subject
                  # Z: Take rows and column corresponding to subject
                  # W: Take rows and columns corresponding to subject
                  
                  sigma.g=sigma.g, Phi=Phi[subday,subday], Psi = matrix(Psi[i+1,i+1]), Sigma =Sigma[idx,idx],
                  X = X[idx,], Z = matrix(Z[idx,(i+1)], ncol = 1), W = W[idx, subday],
                  beta = beta, y = y[idx])
    gamma <- c(gamma, est$par)
    
    # Take Hessian of negative joint likelihood in optimal value of gamma
    H <- c(H, hessian(func = joined.likelihood, x = est$par, 
                      #sigma.g=sigma.g, Phi=Phi[1:days,1:days], 
                      sigma.g=sigma.g, Phi=Phi[subday,subday], 
                      Psi = matrix(Psi[i+1,i+1]), Sigma = Sigma[idx,idx],
                      X = X[idx,], Z = matrix(Z[idx,(i+1)], ncol = 1), W = W[idx, subday],
                      beta = beta, y = y[idx]))
    # H -> negative hessian of joint likelihood
  }
  
  # Get negative joint likelihood using all subject specific gammas
  l.u <- joined.likelihood.H(gamma = gamma, sigma.g = sigma.g, 
                           Phi=Phi, Psi=Psi, Sigma=Sigma, 
                           X=X, Z=Z, W=W, beta=beta, y = y)
  
  #Negative log-likelihood
  #obj <-  l.u + 1/2*logdet(H/(2*pi))
  # Use equation 5.101 to get laplace approximated marginal likelihood of clo
  obj <- l.u + 1/2*sum(log(H/(2*pi)))
  return(obj)
}
opt.fun3(c(0,0,0,0,0,0,0), X = X, Z = Z, W = W, y = y)
par_est3 <- nlminb(c(0,0,0,0,0,0,0), objective = opt.fun3, X = X, Z = Z, W = W, y = y, control = list(trace = 1))
par_est3$objective
# Compare betas
par_est3$par[6:7]

sds <- sqrt(diag(solve(hessian(func = opt.fun3,x=par_est3$par, X = X, Z = Z, W = W, y = y))))

save(par_est, par_est1, par_est2, par_est3, par_est4, file = "all_parameter_estimates.Rdata")

