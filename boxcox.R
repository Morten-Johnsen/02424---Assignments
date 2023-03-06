Trans.eq1 <- function(lambda, y = dioxin$CO){
  y_lambda <- ((y)^lambda - 1)/lambda#, lambda > 0
  return(y_lambda)
}
lambda_NLL <- function(theta, y = dioxin$CO){
  lambda <- theta[1]
  y_lambda <- Trans.eq1(lambda, y)
  mu <- theta[2]#mean(y_lambda)
  sigma <- theta[3]#sd(y_lambda)
  NLL <- -sum(-1/2*log(sigma^2) - (y_lambda-mu)^2 / (2*sigma^2) + (lambda - 1)*log(y))
  return(NLL)
}
theta.hat <- nlminb(start=c(1, 0, 1), objective = lambda_NLL)
lambda.hat <- round(theta.hat$par[1], 3)
mu.hat <- round(theta.hat$par[2], 3)
sd.hat <- round(theta.hat$par[3], 3)

#See if the boxcox transformation is significantly different from the log-transformation
pf_lambda <- function(lambda, y = dioxin$CO){
  n <- length(y)
  y_lambda <- Trans.eq1(lambda)
  sd <- sd(y_lambda)
  ll <- -n/2 * log(sd^2) - n/2 + (lambda - 1)* sum(log(y))
  return(-ll)
}
lambda_interval <- seq(-.5,0.1,0.003)
pf_curve <- -sapply(lambda_interval, FUN = pf_lambda)-max(-sapply(lambda_interval, FUN = pf_lambda))
plot(lambda_interval
     ,pf_curve, type = "l", main = TeX("Profile Likelihood for \\lambda")
     ,ylim = c(-3,0))
grid()
abline(v = lambda.hat, lty = 2)
alpha <- 0.05
c <- -0.5 * qchisq(1-alpha, df = 1)
abline(h = c, col = "red")

max(lambda_interval[pf_curve > c])
lambda.hat
