#### Part 2 ####
library(lme4)
fit0 <- lmer(clo~sex+(1|subjId),data=c.data,REML=FALSE)
summary(fit0)

X <- cbind(1, as.numeric(c.data$sex == "male"))
Z <- dummy(c.data$subjId, levelsToKeep = unique(c.data$subjId))
y <- c.data$clo
dim(Z)
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

#### Mere simpel tilgang med bare at skrive density og likelihood op

opt.fun <- function(theta){
  Psi <- diag(rep(exp(theta[1]),dim(Z)[2]))
  Sigma <- diag(rep(exp(exp(theta[2])),length(y)))
  
  beta <- matrix(theta[3:4], ncol = 1)
  u <- matrix(theta[-c(1:4)], ncol = 1)
  
  V <- Sigma + Z %*% Psi %*% t(Z)
  
  obj <- -1/2*log(det(V))-0.5*t(y - X%*%beta - Z%*%u) %*% solve(V) %*% (y - X%*%beta - Z%*%u)
  
  # fU <- function(Psi, u){
  #   ll <- sum(log(1/(det(Psi) * exp(-0.5*t(u)%*%solve(Psi)%*%u))))
  #   return(-ll)
  # }
  # 
  # fYu <- function(beta, Sigma, u, X, Z, y){
  #   ll <- sum(log(1/det(Sigma) * exp(t(y-X%*%beta-Z%*%u)%*%solve(Sigma)%*%(y - X%*%beta - Z%*%u))))
  # }
  
  #obj <- -(fYu(beta, Sigma, u, X, Z, y) + fU(Psi, u))
  return(obj)
}

opt.fun(theta = c(log(0.11665), log(0.09883), 0.59176, -0.08322, as.vector(t(ranef(fit0)[['subjId']]))))
