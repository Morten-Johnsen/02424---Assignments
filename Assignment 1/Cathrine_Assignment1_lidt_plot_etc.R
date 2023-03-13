# http://web.vu.lt/mif/a.buteikis/wp-content/uploads/PE_Book/3-2-OLS.html

# y = DIOX/log_DIOX
# De variable der giver mening at have med

fit2 <- lm(DIOX_boxcox ~ PLANT + TIME + LAB + O2 + NEFFEKT + QRAT, data = dioxin)
summary(fit2)

dioxin$predicted <- predict(fit2)   # Save the predicted values
dioxin$residuals <- residuals(fit2) # Save the residual values

dioxin %>% 
  gather(key = "iv", value = "x", -DIOX_boxcox, -predicted, -residuals, -LAB, -PLANT, -PRSEK, -OXYGEN, -LOAD, -logDiox, -O2COR, -LOAD_Ordinal, -TIME, -OBSERV, -OXYGEN_Ordinal, -PRSEK_Ordinal, -PLANT_RENO_N, -PLANT_KARA, -LAB_USA_or_KK) %>%  # Get data into shape
  ggplot(aes(x = x, y = DIOX_boxcox)) +  # Note use of `x` here and next line
  geom_segment(aes(xend = x, yend = predicted), alpha = .2) +
  geom_point(aes(color = residuals)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  guides(color = FALSE) +
  geom_point(aes(y = predicted), shape = 1) +
  facet_grid(~ iv, scales = "free_x") +  # Split panels here by `iv`
  theme_bw()



# Mere plot ---------------------------------------------------------------

# Studentized residual is simply a residual divided by its estimated standard deviation.
library(MASS)

#calculate studentized residuals
stud_resids <- studres(fit2)

# Add studendized residulas back to data
dioxin <- cbind(dioxin, stud_resids)

# Sort in order to detect the observations that could be outliers
dioxin[order(-stud_resids),]

# Plot residuals and QQplot
layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = T))
plot(
  fit2$fitted,
  rstudent(fit2),
  main = "Multi Fit Studentized Residuals",
  xlab = "Predictions",
  ylab = "Studentized Resid",
  ylim = c(-2.5, 2.5)
)
abline(h = 0, lty = 2)
plot(
  dioxin$DIOX_boxcox,
  fit2$residuals,
  main = "Residuals by PLANT + TIME + LAB + O2 + NEFFEKT + QRAT",
  xlab = "Predictions (PLANT + TIME + LAB + O2 + NEFFEKT + QRAT)",
  ylab = "Residuals"
)
abline(h = 0, lty = 2)
hist(fit2$resid, main = "Histogram of Residuals")
qqnorm(fit2$resid)
qqline(fit2$resid)





# 
# 
# 
# y = dioxin$H2O
# x = dioxin$O2 
# 
# N = nrow(dioxin)
# beta_0 = 16
# beta_1 = 0.5
# 
# beta_1_est <- cov(x, y) / var(x)
# beta_0_est <- mean(y) - beta_1_est * mean(x)
# print(paste0("Estimated beta_0 = ", beta_0_est, ". True beta_0 = ", beta_0))
# 
# x_mat <- cbind(1, x)
# beta_mat <- solve(t(x_mat) %*% x_mat) %*% t(x_mat) %*% y
# row.names(beta_mat) <- c("beta_0", "beta_1")
# print(beta_mat)
# 
# y_fit <- beta_mat[1] + beta_mat[2] * x
# 
# 
# plot(x = x, y = y, pch = 19, col = "black",yaxs = "i", #Y axis starts at 0, ends at max(Y)
#      xlab = "X", ylab = "Y") # Label the axis
# title("Scatter diagram of (X,Y) sample data and the regression line")
# lines(x = x, y = y_fit, col = "darkgreen", lty = 2)
# # Add Axis labels and ticks at specific positions:
# axis(side = 1, at = c(x[8], x[12]))
# text(x = c(x[8], x[12]), y = -0.2, pos = 1, xpd = TRUE,
#      labels = c(expression(x[8]), expression(x[12])))
# #
# #
# # Add vertical lines:
# segments(x0 = x[8], y0 = 0, x1 = x[8], y1 = y_fit[8], lty = 2)
# segments(x0 = x[8], y0 = 0, x1 = x[8], y1 = y[8])
# segments(x0 = x[12], y0 = 0, x1 = x[12], y1 = y[12])
# #
# #
# # Add some curly brackets:
# pBrackets::brackets(x1 = x[8], y1 = min(y_fit), x2 = x[8], y2 = y[8],
#                     lwd = 1, h = 0.33, col = "red")
# pBrackets::brackets(x1 = x[8], y1 = y_fit[8], x2 = x[8], y2 = min(y_fit),
#                     lwd = 1, h = 0.33, col = "red")
# pBrackets::brackets(x1 = x[12], y1 = y_fit[12], x2 = x[12], y2 = min(y_fit),
#                     lwd = 1, h = 0.33, col = "red")
# pBrackets::brackets(x1 = x[12], y1 = y[12], x2 = x[12], y2 = y_fit[12],
#                     lwd = 1, h = 0.33, ticks = 0.25, col = "red")
# #
# #
# # Add Actual, Fitted and Residual indicator text:
# text(x = 9.9, y = 10.8, labels = expression(Y[8]))
# text(x = 10.9, y = 10.2 , labels = expression(hat(Y)[8]))
# text(x = 13.5, y = 9.1,
#      labels = expression(hat(Y)[12] == "fitted value"))
# text(x = 13.3, y = y_fit[12] + (y[12] - y_fit[12]) * 0.8,
#      labels = expression(hat(epsilon)[12] == "residual"))
# #
# #
# # Add Regression line
# text(x = x[9], y = y[11],
#      labels = expression(hat(Y) == hat(beta)[0] + hat(beta)[1] * X))
# arrows(x0 = x[9], y0 = y_fit[10]*1.05, x1 = x[10], y1 = y_fit[10], length = 0.1, col = "red")
