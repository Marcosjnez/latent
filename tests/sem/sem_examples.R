library(latent)
library(lavaan)

model <- '
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
set.seed(2026)
estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

fit <- lcfa(model = model,
            data = PoliticalDemocracy,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = FALSE,
            do.fit = TRUE)
fit@modelInfo$param
latInspect(fit, what = "convergence")
fit@Optim$elapsed
fit@transformed_pars$B
fit@modelInfo$control_optimizer$deltaparam
fit@parameters

fit2 <- sem(model, data = PoliticalDemocracy,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood)
# Compare the ML discrepancy and the log likelihoods of the same three models.
reference <- lavaan::fitMeasures(fit2, c("fmin", "logl", "unrestricted.logl",
                                         "baseline.chisq"))
ll <- latent::latInspect(fit, "loglik")[, "overall"]
comparison <- rbind(
  latent = c(latent::latInspect(fit, "loss")["loss", "overall"],
             ll["loglik"], ll["loglik_base"], ll["loglik_sat"]),
  lavaan = c(2*reference["fmin"], reference["logl"],
             reference["unrestricted.logl"] - reference["baseline.chisq"]/2,
             reference["unrestricted.logl"]))
colnames(comparison) <- c("loss", "loglik", "loglik_base", "loglik_sat")
comparison
stopifnot(max(abs(comparison[1L, ] - comparison[2L, ])) < 1e-5)

lavaan::lavInspect(fit2, "est")$beta
round(latent::latInspect(fit, "est", sort = FALSE)$B, 3)
fit2
