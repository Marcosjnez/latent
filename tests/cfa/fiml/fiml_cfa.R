library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

data_missing <- HolzingerSwineford1939
data_missing[1:5, "x1"] <- NA
data_missing[9:20, "x5"] <- NA
data_missing[167:189, "x6"] <- NA
data_missing[244:251, "x9"] <- NA
data_missing[280, "x8"] <- NA
# rvars <- sample(paste("x", c(1, 4, 6), sep = ""), size = nmiss, replace = TRUE)
# rsubjs <- sample(1:1000, size = nmiss)
# for(j in 1:nmiss) data_missing[rsubjs[j], rvars[j]] <- NA

estimator <- "ml"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

set.seed(2026)
fit <- lcfa(data = data_missing,
            model = model,
            estimator = estimator,
            missing = "ml",
            group = "school",
            meanstructure = meanstructure,
            std.lv = std.lv,
            std.ov = std.ov,
            likelihood = likelihood,
            # missing = "pairwise.complete.obs",
            positive = FALSE,
            se = "standard",
            do.fit = TRUE)

# fit@modelInfo$param

# With lavaan:
fit2 <- cfa(model = model,
            data = data_missing,
            estimator = estimator,
            missing = "fiml",
            group = "school",
            meanstructure = meanstructure,
            std.lv = std.lv,
            std.ov = std.ov,
            likelihood = likelihood)

# The direct pattern likelihood has a different additive constant in its
# optimization loss. Compare the observed-data log likelihood on the same cases.
ll <- latInspect(fit, "loglik")["loglik", "overall"]
c(latent = ll, lavaan = as.numeric(logLik(fit2)))
stopifnot(abs(ll - as.numeric(logLik(fit2))) < 1e-5)

latInspect(fit, "est", sort = FALSE)
lavInspect(fit2, "est")

lavInspect(fit2, "se")
fit2@ParTable$se
blocks <- with(fit@dataList$data_param, c(lambda_group, theta_group, psi_group))
x <- se(fit, parameters = fit@modelInfo$trans[blocks])
x$table_se
