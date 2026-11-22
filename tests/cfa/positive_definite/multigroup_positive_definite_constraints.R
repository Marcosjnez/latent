library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9
          x1 ~~ x5
          x1 ~~ x4
          x4 ~~ x5
          x4 ~~ x6'

# With latent:
estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

# Obtain starting values from the unrestricted covariance fit. lavaan
# warns about a negative residual variance in this example.
fit2 <- lavaan::cfa(model, data = HolzingerSwineford1939,
                    group = "school", estimator = estimator,
                    std.lv = std.lv, std.ov = std.ov,
                    meanstructure = meanstructure, likelihood = likelihood)
est <- lavaan::lavInspect(fit2, "est")
start <- list()
for(g in names(est)) {
  theta <- est[[g]]$theta
  smallest <- min(eigen(theta, symmetric = TRUE, only.values = TRUE)$values)
  diag(theta) <- diag(theta) + max(0, 0.1 - smallest)
  start[[paste0("xtheta.", g)]] <- chol(theta)
  start[[paste0("xpsi.", g)]] <- chol(est[[g]]$psi)
  start[[paste0("lambda.", g)]] <- est[[g]]$lambda
}

set.seed(2026)
fit <- lcfa(data = HolzingerSwineford1939,
            model = model,
            group = "school",
            estimator = estimator,
            ordered = FALSE,
            std.ov = std.ov,
            std.lv = std.lv,
            positive = TRUE,
            penalties = TRUE,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = "standard",
            do.fit = TRUE,
            control = list(opt = "newton", step = "trust", tcg_maxit = 100L,
                           start = start, rstarts = 1L, se_method = "RH"))
stopifnot(fit@Optim$convergence)
# Expected log likelihood is about -3674.300. Redundant Gram coordinates
# can produce variance warnings even at a converged fit.
latInspect(fit, "loglik")

fit@Optim$iterations
fit@Optim$convergence
fit@timing

# Compare the unconstrained fit:
fitmeasures(fit2, fit.measures = c("cfi", "tli", "rmsea", "srmr"))
fitMeasures(fit2, "logl")
fitMeasures(fit2, "unrestricted.logl")
fitMeasures(fit2, "baseline.chisq")

diag(lavaan::lavInspect(fit2, "est")$Pasteur$theta) # Unconstrained residual variances
det(inspect(fit2, what = "est")$Pasteur$theta)
