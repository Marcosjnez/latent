library(latent)
gss <- latent::gss82
gss$EDUCR <- as.integer(gss$EDUCR)-1L
set.seed(20216)
fit <- lca(data = gss,
           nclasses = 3L,
           multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
           # covariates = c("RACE", "SEX", "EDUCR", "AGE"),
           # outcomes = "MARITAL",
           # adjustment = "bk",
           # classification = "modal",
           # model = list("UNDERSTA ~~ COOPERAT
           #               PURPOSE ~~ COOPERAT"),
           # start = start,
           penalties = list(class = list(alpha=1),
                            prob  = list(alpha=0)),
           # control = list(opt = "lbfgs", step = "armijo",
           #                rstarts = 30, cores = 30),
           control = list(opt = "em", maxit = 1000L, eps = 1e-10,
                          mopt = "lbfgs", mstep_maxit = 50L, mstep_eps = 1e-05,
                          step = "armijo", step_maxit = 30L,
                          rstarts = 30, cores = 30),
           do.fit = TRUE)
latInspect(fit, what = "loglik")
latInspect(fit, what = "convergence")
latInspect(fit, what = "elapsed")
# loglik: -3891.252 # penalized_loglik: -3892.478
# With both residual dependencies and direct L-BFGS:
# loglik: -3879.167 # penalized_loglik: -3880.371

jacob <- jacobian(fit, parameters = fit@modelInfo$trans[names(fit@modelInfo$param)])
dim(jacob)
labs1 <- c(fit@modelInfo$trans$PURPOSE)
labs2 <- c(fit@modelInfo$trans$PURPOSE)
jacob[labs1, labs2]

dconstr <- constraints_derivs(fit, parameters = fit@modelInfo$trans[names(fit@modelInfo$param)])
dim(dconstr$dconstr)
labs1 <- c(fit@modelInfo$trans$ACCURACY)
dconstr$dconstr[labs1, , drop = FALSE]

# Print model fit info:
fit

# Get fit indices:
getfit(fit)

# Get a summary:
summary(fit)

# Bivariate residuals:
lbvr(fit)

# Plot:
plot(fit)
plot(fit, type = "coefficients", what = "OR", se_type = "information")

# With the covariates above, predict class membership for two respondents:
if(length(fit@dataList$covariates_names)) {
  df <- data.frame(RACE = c("BLACK", "WHITE"),
                   SEX = c("MALE", "MALE"),
                   EDUCR = c(3, 3),
                   AGE = c(35, 35))
  predict(fit, new = df)
} else {
  head(predict(fit)) # Common class probabilities in the intercept-only model
}

# Inspect model objects:
latInspect(fit, what = "convergence")
latInspect(fit, what = "profile")
latInspect(fit, what = "coefs")
latInspect(fit, what = "classification")
latInspect(fit, what = "pattern")
latInspect(fit, what = "table")

# Get standard errors:
SE <- se(fit, type = "information", digits = 4,
         parameters = fit@modelInfo$trans[names(fit@modelInfo$param)])
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "information", confidence = 0.95, digits = 2)
CI$table
