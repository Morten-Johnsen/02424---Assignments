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
fit0 <- lmer(clo~sex+(1|subjId),data=c.data,REML=FALSE)
summary(fit0)

X <- cbind(1, as.numeric(c.data$sex == "male"))
Z <- dummy(c.data$subjId, levelsToKeep = unique(c.data$subjId))
y <- c.data$clo
dim(Z)

#### Mere simpel tilgang med bare at skrive density og likelihood op
## Joint Likelihood
library(mvtnorm)
library(numDeriv)
inner.nll <- function(u,beta,sigma,sigma.u,X,Z,y){
  u <- matrix(u, ncol = 1)
  return(-sum(dnorm(y, mean = X%*%beta + Z%*%u, sd = sqrt(sigma), log = T), 
              dnorm(u, sd = sqrt(sigma.u), log = T)))
  # return(-dmvnorm(y, mean = X%*%beta + Z%*%u, sigma = Sigma, log = T)
  #        -dmvnorm(t(u), sigma = Psi, log = T))
}

## inner-outer
grad <- function(u, beta, X, Z, sigma.u, sigma, y){
  Psi <- diag(rep(sigma.u,dim(Z)[2]))
  Sigma <- diag(rep(sigma,length(y)))
  return(t(Z)%*%solve(Sigma)%*%(y - X%*%beta - Z%*%u)-solve(Psi)%*%u)
}

hess <- function(u, beta, X, Z, sigma.u, sigma, y){
  Psi <- diag(rep(sigma.u,dim(Z)[2]))
  Sigma <- diag(rep(sigma,length(y)))
  return(-t(Z)%*%solve(Sigma)%*%Z-solve(Psi))
}

nll <- function(theta, X, Z, y, save_u = F){
  beta <- matrix(theta[1:2], ncol = 1)
  sigma.u <- exp(theta[3])
  sigma <- exp(theta[4])
  
  est <- nlminb(rep(0, dim(Z)[2]),objective = inner.nll, #gradient = grad, hessian = hess,
                beta=beta, sigma.u=sigma.u, sigma=sigma, X=X, Z = Z, y = y)
  print(est$objective)
  u <- est$par
  if (save_u){
    u <<- est$par
  }
  l.u <- est$objective
  
  #H <- diag(hessian(func = inner.nll, x = est$par, beta = beta, sigma = sigma, sigma.u = sigma.u, X = X, Z = Z, y = y))
  H <- diag(hess(u, beta, X, Z, sigma.u, sigma, y))
  return(l.u + 0.5 * sum(log(abs(H/(2*pi)))))
}
#library(profvis)
#profvis(nll(c(0.59176, -0.08322, log(0.1),log(0.13)), X, Z, y))
par_est_inner_opt <- nlminb(start = c(0, 0, 0, 0), objective = nll, X = X, Z = Z, y = y, control = list(trace=1))

exp(par_est_inner_opt$par[3:4])
nll(par_est_inner_opt$par, X, Z, y, save_u = T)
cbind(u, as.numeric(t(ranef(fit0)$subjId)))

opt.fun <- function(theta){
  Psi <- diag(rep(exp(theta[1]),dim(Z)[2]))
  Sigma <- diag(rep(exp(theta[2]),length(y)))
  
  beta <- matrix(theta[3:4], ncol = 1)
  
  V <- Sigma + Z%*%Psi%*%t(Z)
  
  #log-likelihood of multivariate normal distribution:
  obj <-  mvtnorm::dmvnorm(y, mean = X%*%beta, sigma = V, log = T)
  
  return(-obj)
}

par_est <- nlminb(start = c(log(0.11665),log(0.09883), 0.59176, -0.08322),
                  objective = opt.fun, control = list(trace = 1))
sqrt(exp(par_est$par[1:2]))
sds <- sqrt(diag(solve(hessian(func = opt.fun, x = par_est$par))))

manual_fit0 <- data.frame("Parameter" = c("Psi (Sigma.u)", "Sigma", "beta0", "beta1"),
                          "Interpretation" = c("SubjId", "Model Residual", "Model intercept", "Model slope (sex)"),
                          "Estimate" = c(exp(par_est$par[1:2]), par_est$par[3:4]))
cat("Log-likelihood = ", -par_est$objective)
manual_fit0
summary(fit0)

Psi <- diag(rep(exp(par_est$par[1]),dim(Z)[2]))
Sigma <- diag(rep(exp(par_est$par[2]),length(y)))
beta <- matrix(par_est$par[3:4], ncol = 1)

u_est <- function(u, beta, Psi, Sigma, X, Z, y){
  u <- matrix(u, ncol=1)
  obj <- -1/2*log(max(det(Sigma), 1e-300)) -1/2*t(y- X%*%beta-Z%*%u)%*%solve(Sigma)%*%(y - X%*%beta - Z%*%u) - 1/2*log(max(det(Psi), 1e-300)) - 1/2*t(u)%*%solve(Psi)%*%u
  return(-obj)
}

u_par <- nlminb(start = rep(0,47), objective = u_est, beta = beta, Psi = Psi, Sigma = Sigma, 
                X = X, Z = Z, y = y, control = list(trace = 1))

y_adj <- y - X%*%beta
u <- solve(t(Z)%*%solve(Sigma)%*%Z + solve(Psi)) %*% (t(Z)%*%solve(Sigma)%*%y_adj)
cbind(ranef(fit0)$subjId, u, u_inner, u_par$par)

#### 2.2 ####
fit1 <- lmer(clo~sex+(1|subjId)+(1|subjId:day),data=c.data,REML=FALSE)
ranef(fit1)

X <- cbind(1, as.numeric(c.data$sex == "male"))
Z <- dummy(c.data$subjId, levelsToKeep = unique(c.data$subjId))
#Introduce W which is the design matrix for the nested subjdID:day observations
c.data$subjDay <- factor(paste0(c.data$subjId, "/",c.data$day))
W <- dummy(c.data$subjDay, levelsToKeep = unique(c.data$subjDay))
y <- c.data$clo
dim(Z)

#### Mere simpel tilgang med bare at skrive density og likelihood op

opt.fun1 <- function(theta){
  Psi <- diag(rep(exp(theta[2]),dim(Z)[2]))
  Sigma <- diag(rep(exp(theta[3]),length(y)))
  Phi <- diag(rep(exp(theta[1]),dim(W)[2]))
  
  beta <- matrix(theta[4:5], ncol = 1)
  
  V <- Sigma + Z%*%Psi%*%t(Z) + W%*%Phi%*%t(W)
  
  obj <-  mvtnorm::dmvnorm(y, mean = X%*%beta, sigma = V, log = T)
  
  return(-obj)
}

par_est1 <- nlminb(start = c(0, 0, 0, 1, 1),
                   objective = opt.fun1, control = list(trace = 1))
#Sammenlign med fit1
sds1 <- sqrt(diag(solve(hessian(func = opt.fun1, x = par_est1$par))))
summary(fit1)
manual_fit1 <- data.frame("Parameter" = c("Phi (Sigma.v)", "Psi (Sigma.u)", "Sigma", "beta0", "beta1"),
                          "Interpretation" = c("SubjId/day", "SubjId", "Model Residual", "Model Intercept", "Model Slope (sex)"),
                          "Estimate" = c(exp(par_est1$par[1:3]), par_est1$par[4:5]))
cat("Log-likelihood = ", -par_est1$objective)
manual_fit1
sds1

#Hierarchical estimation: We first estimated Beta, Psi and Phi.
#now we will estimate u and the nested value v.
Phi <- diag(rep(exp(par_est1$par[1]),dim(W)[2]))
Psi <- diag(rep(exp(par_est1$par[2]),dim(Z)[2]))
Sigma <- diag(rep(exp(par_est1$par[3]),length(y)))
beta <- matrix(par_est1$par[4:5], ncol = 1)

A <- rbind(
  cbind(t(Z)%*%solve(Sigma)%*%Z+solve(Psi), t(Z)%*%solve(Sigma)%*%W),
  cbind(t(W)%*%solve(Sigma)%*%Z, t(W)%*%solve(Sigma)%*%W+solve(Phi))
)
B <- rbind(t(Z)%*%solve(Sigma)%*%(y-X%*%beta),
           t(W)%*%solve(Sigma)%*%(y-X%*%beta))
uv <- solve(A, B)

u <- t(t(uv[1:dim(Z)[2]]))
v <- t(t(uv[-c(1:dim(Z)[2]),]))

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

cbind(ranef(fit1)$subjId, u)
cbind(ranef(fit1)$`subjId:day`, v)

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

idx_Psi <- subjSex$sex
idx_Phi <- subjDaySex$sex
idx_Sigma <- as.numeric(c.data$sex == "male")


opt.fun2 <- function(theta){
  alpha <- theta[1] #constant effect of gender across all variances
  
  #Estimation in the exponential domain and estimate alpha as a scaling factor of the female variance
  # -> e(theta)*e(alpha) = e(theta + alpha)
  Psi <- diag(exp(theta[3] + idx_Psi*alpha))
  Sigma <- diag(exp(theta[4] + idx_Sigma*alpha))
  Phi <- diag(exp(theta[2] + idx_Phi*alpha))
  
  beta <- matrix(theta[5:6], ncol = 1)
  
  V <- Sigma + Z%*%Psi%*%t(Z) + W%*%Phi%*%t(W)
  
  obj <-  mvtnorm::dmvnorm(y, mean = X%*%beta, sigma = V, log = T)
  
  return(-obj)
}

par_est2 <- nlminb(start = c(0, 0, 0, 0, 1, 1),
                   objective = opt.fun2, control = list(trace = 1))
sds2 <- sqrt(diag(solve(hessian(opt.fun2, x = par_est2$par))))
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


#### 2.5 ####

#### 2.6 ####
#For this expression we do not know the marginal likelihood and instead have to approximate it by
#using laplace approximation.
#The inner loop is the estimation of u,v and gamma and should therefore use the joint likelihood.

grad2 <- function(uv, beta, sigma, sigma.u, sigma.v, X, Z, W, y){
  #for diagonal matrices: inverse(M) = 1/M
  u_length <- dim(Z)[2]
  u <- matrix(uv[1:u_length], ncol = 1)
  v <- matrix(uv[-c(1:u_length)], ncol = 1)
  invPsi <-   diag(1/sigma.u, nrow = u_length, ncol = u_length)
  invPhi <-   diag(1/sigma.v, nrow = dim(W)[2], ncol = dim(W)[2])
  invSigma <- diag(1/sigma, nrow = length(y), ncol = length(y))
  
  e <- y - X%*%beta - Z%*%u - W%*%v
  invSigma_x_e <- invSigma%*%e
  gradient <- matrix(c(t(Z)%*%invSigma_x_e-invPsi%*%u, t(W)%*%invSigma_x_e-invPhi%*%v), ncol = 1)
  return(-gradient)
}

hess2 <- function(uv, beta, sigma, sigma.u, sigma.v, X, Z, W, y){
  #for diagonal matrices: inverse(M) = 1/M
  u_length <- dim(Z)[2]
  u <- matrix(uv[1:u_length], ncol = 1)
  v <- matrix(uv[-c(1:u_length)], ncol = 1)
  invPsi <-   diag(1/sigma.u, nrow = u_length, ncol = u_length)
  invPhi <-   diag(1/sigma.v, nrow = dim(W)[2], ncol = dim(W)[2])
  invSigma <- diag(1/sigma, nrow = length(y), ncol = length(y))
  
  #Save one matrix multiplication pr. run
  tWinvSigma <- t(W)%*%invSigma
  #Save one more matrix multiplication pr. run
  tZinSigma <- t(Z)%*%invSigma
  
  luu <- -tZinSigma%*%Z-invPsi
  lvu <- -tWinvSigma%*%Z
  luv <- -tZinSigma%*%W
  lvv <- -tWinvSigma%*%W-invPhi
  
  hessian <- rbind(cbind(luu, luv),
                   cbind(lvu, lvv))
  return(-hessian)
}

library(msos)
joint.likelihood <- function(uv, beta, sigma, sigma.u, sigma.v, X, Z, W, y){
  u_length <- dim(Z)[2]
  u <- matrix(uv[1:u_length], ncol = 1)
  v <- matrix(uv[-c(1:u_length)], ncol = 1)

  jll <- sum(dnorm(y, mean = X%*%beta + Z%*%u + W%*%v, sd = sqrt(sigma), log = T)) +
                           sum(dnorm(as.numeric(u), sd = sqrt(sigma.u), log = T)) + 
                           sum(dnorm(as.numeric(v), sd = sqrt(sigma.v), log = T))
  return(-jll)
}

nll2 <- function(theta, X, Z, W, y, save_u = F){
  beta <- matrix(theta[1:2], ncol = 1)
  sigma.u <- exp(theta[3])
  sigma.v <- exp(theta[4])
  sigma <- exp(theta[5])
  
  Psi <-   diag(rep(sigma.u,dim(Z)[2]))
  Phi <-   diag(rep(sigma.v,dim(W)[2]))
  Sigma <- diag(rep(sigma,length(y)))
  
  est <- nlminb(start = rep(0, dim(Z)[2]+dim(W)[2]), objective = joint.likelihood,
                beta = beta, sigma = sigma, sigma.u = sigma.u, sigma.v = sigma.v,
                X =X, Z = Z, W = W, y = y)
                # gradient = grad2,
                # hessian = hess2)
  uv <- est$par
  l.u <- est$objective
  #Faster solution: Find the optimal solution based on solving a system of equations (have to know the hessian)
  # A2 <- rbind(
  #   cbind(t(Z)%*%solve(Sigma)%*%Z+solve(Psi), t(Z)%*%solve(Sigma)%*%W),
  #   cbind(t(W)%*%solve(Sigma)%*%Z, t(W)%*%solve(Sigma)%*%W+solve(Phi))
  # )
  # B2 <- rbind(t(Z)%*%solve(Sigma)%*%(y-X%*%beta),
  #             t(W)%*%solve(Sigma)%*%(y-X%*%beta))
  # uv2 <- solve(A2, B2)
  # u <- uv2[1:dim(Z)[2]]
  # v <- uv2[-c(1:dim(Z)[2])]
  
  #evaluate the joint likelihood of this solution
  # l.u <- -(mvtnorm::dmvnorm(y, mean = X%*%beta + Z%*%u + W%*%v, sigma = Sigma, log = T) +
  #            mvtnorm::dmvnorm(as.numeric(u), sigma = Psi, log = T) +
  #            mvtnorm::dmvnorm(as.numeric(v), sigma = Phi, log = T))
  
  #Return the NEGATIVE hessian according to the section between 5.100 and 5.101 p. 201
  H <- -hess2(uv, beta = beta, sigma = sigma, sigma.u = sigma.u, sigma.v = sigma.v, 
              X = X, Z = Z, W = W, y = y)
  #return negative log likelihood (and thus + 0.5log(prod....))
  return(l.u + 0.5*logdet(H/(2*pi)))
}
library(profvis)
profvis(nll2(theta = c(0.59242, -0.08439, log(.09730^2),log(.10395^2), log(.05597^2)),
             X = X, Z = Z, W = W, y = y))
par_est_inner_opt2 <- nlminb(start = c(0.59242, -0.08439, log(.09730^2),log(.10395^2), log(.05597^2)), 
                             objective = nll2, 
                             X = X, Z = Z, W = W, y = y, 
                             control = list(trace=1))

sqrt(exp(par_est_inner_opt2$par[3:5]))
exp(par_est2$par[1:3])

hessian(func = joint.likelihood, x = c(0.59242, -0.08439, log(.09730^2),log(.10395^2), log(.05597^2)),
                                beta, sigma, sigma.u, sigma.v, X, Z, W, y)
        