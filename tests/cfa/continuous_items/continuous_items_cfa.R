library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'
S <- cov(HolzingerSwineford1939[, paste("x", 1:9, sep = "")])
means <- colMeans(HolzingerSwineford1939[, paste("x", 1:9, sep = "")])

set.seed(2026)
estimator <- "uls"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

fit <- lcfa(model = model,
            data = HolzingerSwineford1939,
            # sample.cov = S,
            # sample.mean = means,
            # sample.nobs = 301,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = "standard",
            # control = list(opt = "lbfgs", step = "armijo"),
            do.fit = TRUE)
latInspect(fit, what = "convergence")
fit@Optim$elapsed
fit@modelInfo$param

# fit@modelInfo$param

# With lavaan:
fit2 <- lavaan::cfa(model = model,
                    data = HolzingerSwineford1939,
                    # sample.cov = S,
                    # sample.mean = means,
                    # sample.nobs = 301,
                    estimator = estimator,
                    std.lv = std.lv,
                    std.ov = std.ov,
                    meanstructure = meanstructure,
                    likelihood = likelihood,
                    do.fit = TRUE)
# The default ULS criteria use different covariance divisors and weights.
# latent uses S * (n-1)/n and counts both triangles in its squared residuals.
latent_loss <- latInspect(fit, "loss")["loss", "overall"]
c(latent = latent_loss, lavaan_twice_fmin = 2 * unname(fitMeasures(fit2, "fmin")))

# Supply the same sample moments and weights for a direct comparison.
# lavaan stacks the means followed by the unique entries of the covariance.
n <- nrow(HolzingerSwineford1939)
W <- matrix(1, nrow(S), ncol(S))
diag(W) <- 0.5
WLS.V <- diag(c(rep(0.5, length(means)), W[lower.tri(W, diag = TRUE)]))
fit2_matched <- lavaan::cfa(model, sample.cov = S * (n-1)/n, sample.mean = means,
                            sample.nobs = n, sample.cov.rescale = FALSE,
                            ordered = FALSE, estimator = "DWLS", WLS.V = WLS.V,
                            std.lv = std.lv, std.ov = std.ov,
                            meanstructure = meanstructure, se = "none")
# lavaan's least-squares test statistic uses n-1 observations.
lavaan_loss <- lavaan::lavInspect(fit2_matched, "test")$standard$stat/(n-1)
c(latent = latent_loss, lavaan_matched = lavaan_loss)
stopifnot(abs(latent_loss - lavaan_loss) < 1e-7)

lavaan::lavInspect(fit2, "est")
latInspect(fit, "est", sort = FALSE)
fit@Optim$SE$se
lavaan::lavInspect(fit2, "se")
fitMeasures(fit2, c("baseline.chisq", "rmr"))
getfit(fit)
