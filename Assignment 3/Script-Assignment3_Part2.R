######################## Part 2 ########################
library(lme4)
#### Part 2.1 ####
fit0 <- lmer(clo~sex+(1|subjId),data=c.data,REML=FALSE)
summary(fit0)

X <- cbind(1, as.numeric(c.data$sex == "male"))
Z <- dummy(c.data$subjId, levelsToKeep = unique(c.data$subjId))
y <- c.data$clo
dim(Z)

#### Mere simpel tilgang med bare at skrive density og likelihood op

opt.fun <- function(theta){
  Psi <- diag(rep(exp(theta[1]),dim(Z)[2]))
  Sigma <- diag(rep(exp(theta[2]),length(y)))
  
  beta <- matrix(theta[3:4], ncol = 1)
  
  V <- Sigma + Z%*%Psi%*%t(Z)
  
  obj <-  mvtnorm::dmvnorm(y, mean = X%*%beta, sigma = V, log = T)
  
  return(-obj)
}

par_est <- nlminb(start = c(log(0.11665),log(0.09883), 0.59176, -0.08322),
                  objective = opt.fun, control = list(trace = 1))
sqrt(exp(par_est$par[1:2]))
manual_fit0 <- data.frame("Parameter" = c("Psi (Sigma.u)", "Sigma", "beta0", "beta1"),
                          "Interpretation" = c("SubjId", "Model Residual", "Model intercept", "Model slope (sex)"),
                          "Estimate" = c(exp(par_est$par[1:2]), par_est$par[3:4]))
cat("Log-likelihood = ", -par_est$objective)
manual_fit0
summary(fit0)

Psi <- diag(rep(exp(par_est$par[1]),dim(Z)[2]))
Sigma <- diag(rep(exp(par_est$par[2]),length(y)))
beta <- matrix(par_est$par[3:4], ncol = 1)

y_adj <- y - X%*%beta
u <- solve(t(Z)%*%solve(Sigma)%*%Z + solve(Psi)) %*% (t(Z)%*%solve(Sigma)%*%y_adj)
cbind(ranef(fit0)$subjId, u)

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
summary(fit1)
manual_fit1 <- data.frame("Parameter" = c("Phi (Sigma.v)", "Psi (Sigma.u)", "Sigma", "beta0", "beta1"),
                          "Interpretation" = c("SubjId/day", "SubjId", "Model Residual", "Model Intercept", "Model Slope (sex)"),
                          "Estimate" = c(exp(par_est1$par[1:3]), par_est1$par[4:5]))
cat("Log-likelihood = ", -par_est1$objective)
manual_fit1

#Hierarchical estimation: We first estimated Beta, Psi and Phi.
#now we will estimate u and the nested value v.
Phi <- diag(rep(exp(par_est1$par[1]),dim(W)[2]))
Psi <- diag(rep(exp(par_est1$par[2]),dim(Z)[2]))
Sigma <- diag(rep(exp(par_est1$par[3]),length(y)))
beta <- matrix(par_est1$par[4:5], ncol = 1)

v <- matrix(0, nrow = dim(W)[2])
u <- matrix(0, nrow = dim(Z)[2])

v_old <- v
u_old <- u
iteration <- 0

while(iteration == 0 | max(abs(v_old - v)) + max(abs(u_old - u)) > 1e-5){
  iteration <- iteration + 1
  u_old <- u
  v_old <- v
  
  u <- solve(t(Z)%*%solve(Sigma)%*%Z + solve(Psi), t(Z)%*%solve(Sigma)%*%(y - X%*%beta - W%*%v))
  v <- solve(t(W)%*%solve(Sigma)%*%W + solve(Phi), t(W)%*%solve(Sigma)%*%(y - X%*%beta - Z%*%u))
  print(iteration)
}

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

manual_fit2 <- data.frame("Parameter" = c("alpha", "Phi (Sigma.v)", "Psi (Sigma.u)", "Sigma", "beta0", "beta1"),
                          "Interpretation" = c("Male gender variance scaling", "SubjId/day", "SubjId", "Model Residual", "Model Intercept", "Model Slope (sex)"),
                          "Estimate" = c(exp(par_est2$par[1:4]), par_est2$par[5:6]))
cat("Log-likelihood = ", -par_est2$objective)
manual_fit2

#### 2.4 ####

#### 2.5 ####

#### 2.6 ####



###### GAMMELT #####
#Simultaneous estimation of beta and u for known variances p. 184
beta <- solve(t(X)%*%X)%*%t(X)%*%y
beta_old <- beta
u <- matrix(0, nrow = 47)
u_old <- u

Sigma <- diag(rep(1,length(y))); Psi <- diag(rep(1,dim(Z)[2]))

iterations <- 0
#Ved ikke om det her er rigtigt - det er som om der mangler et eller andet eller at det kan gøres smartere... 
#Det virker dog til at få de rigtige parameter estimater
while ((all(abs(beta - beta_old) > 1e-9) & all(abs(u - u_old) > 1e-9)) | iterations < 1){
  
  beta_old <- beta
  u_old <- u
  iterations <- iterations + 1
  #calculate the adjusted observation
  y_adj <- y - X%*%beta
  #estimate u
  
  u <- solve(t(Z)%*%solve(Sigma)%*%Z + solve(Psi)) %*% (t(Z)%*%solve(Sigma)%*%y_adj)
  
  y_adj <- y - Z%*%u
  
  #reestimate beta
  beta <- solve(t(X)%*%solve(Sigma)%*%X)%*%t(X)%*%solve(Sigma)%*%y_adj
  
  #tmp.func <- function(theta, u, e, Z, X){
  tmp.func <- function(theta, u, e){
    Sigma <- exp(theta[1])
    Psi <- exp(theta[2])
    
    obj <- sum(dnorm(e, sd = sqrt(Sigma), log = T)) + sum(dnorm(u, sd = sqrt(Psi), log = T))
    #Sigma <- diag(rep(exp(theta[1]),length(y))); Psi <- diag(rep(exp(theta[2]),dim(Z)[2]))
    
    #V <- Sigma + Z %*% Psi %*% t(Z)
    #obj <- -0.5 * log(det(V)) - 0.5*log(det(t(X)%*%solve(V)%*%X)) - 0.5*t(e) %*% solve(V) %*% e
    return(-obj)
  }
  e <- y - X%*%beta - Z%*%u
  est <- nlminb(start = c(1,1), objective = tmp.func, u = u, e = e)
  #est <- nlminb(start = c(-2,-1), objective = tmp.func, u = u, e = e, Z = Z, X = X)
  
  Sigma <- diag(rep(exp(est$par[1]),length(y)))
  Psi <- diag(rep(exp(est$par[2]),dim(Z)[2]))
  
  
  if(iterations %in% c(1, 10, 100, 200, 300, 400, 600, 800, 1000, 10000)){
    cat("\nIteration: ", iterations, " done. Update difference: ",max(abs(beta-beta_old))," \n----------------------------")
  }
}
iterations
beta-beta_old
fit0
beta
cbind(ranef(fit0)$subjId, u)
sqrt(unique(diag(Psi)))
sqrt(unique(diag(Sigma)))
fit0