library(latent)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

fit <- lcfa(HolzingerSwineford1939,
            model = model,
            group = "school",
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            likelihood = likelihood,
            se = "standard",
            control = NULL,
            do.fit = TRUE)

latInspect(fit, "est")
# fit@modelInfo$param

# With lavaan:
fit2 <- lavaan::cfa(data = HolzingerSwineford1939,
                    model = model,
                    group = "school",
                    estimator = estimator,
                    std.ov = std.ov,
                    std.lv = std.lv,
                    meanstructure = meanstructure,
                    likelihood = likelihood,
                    do.fit = TRUE)
loss <- latInspect(fit, "loss")["loss", "overall"]
ll <- latInspect(fit, "loglik")
c(latent = loss, lavaan = 2 * unname(fitMeasures(fit2, "fmin")))
c(latent = ll["loglik", "overall"], lavaan = as.numeric(logLik(fit2)))
stopifnot(abs(loss - 2 * fitMeasures(fit2, "fmin")) < 1e-7,
          abs(ll["loglik", "overall"] - as.numeric(logLik(fit2))) < 1e-5,
          abs(ll["loglik_sat", "overall"] - fitMeasures(fit2, "unrestricted.logl")) < 1e-5)
