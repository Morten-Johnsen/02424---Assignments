rm(list=ls())
library(GGally)
source("DataPrep.R")

Trans.eq1 <- function(lambda, y = dioxin$DIOX){
  y_lambda <- ((y)^lambda - 1)/lambda#, lambda > 0
  return(y_lambda)
}
lambda_NLL <- function(theta, y = dioxin$DIOX){
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
pf_lambda <- function(lambda, y = dioxin$DIOX){
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

#95% profilelikelihood CI for lambda include 0 and thus a log-transformation is sufficient and to be preferred over the boxcox.
#wald CI for confirmation:
sd_reg <- sqrt(diag(solve(hessian(func = lambda_NLL, x = theta.hat$par))))
lambda.hat + qt(c(alpha/2, 1-alpha/2), df = length(dioxin$DIOX) - 1) * sd_reg[1]


dioxin$DIOX_boxcox <- Trans.eq1(lambda.hat, dioxin$DIOX)
hist(dioxin$DIOX_boxcox, breaks = 5)

dioxin %>%
  mutate(logDiox = log(DIOX)) %>%
  select(-PLANT, -LAB, -OXYGEN, -LOAD, -PRSEK, -OBSERV) %>%
  melt() %>%
  ggplot()+
  geom_histogram(aes(x = value), bins = 10)+
  facet_wrap(~variable, scales = "free")

#Block effects: PLANT (3 plants, RENO_N, RENO_S and KARA), TIME (For RENO_N the experiment
#was repeated at a later time point, 2, as well.), LAB (Two labs. One in DK and one in USE)
#considerable measurement noise is expected.

dioxin %>%
  mutate(logDiox = log(DIOX)) %>%
  select(logDiox, TIME, LAB_USA_or_KK
         , PLANT_RENO_N, PLANT_RENO_S, PLANT_KARA
         , OXYGEN_Ordinal
         , LOAD_Ordinal
         , PRSEK_Ordinal
         , O2, O2COR, NEFFEKT, QRAT) %>%
  ggpairs()

#Active variables: 
#   "Theoretical": OXYGEN, LOAD, PRSEK.
#   "Measured": O2 (O2COR), NEFFEKT, QRAT

#### 2) ####
fit <- lm(DIOX_boxcox ~ , data = dioxin)
summary(fit)


#### 3) ####
fit2 <- lm(DIOX_boxcox ~ PLANT + TIME + LAB + O2 + NEFFEKT + QRAT, data = dioxin)
summary(fit2)
