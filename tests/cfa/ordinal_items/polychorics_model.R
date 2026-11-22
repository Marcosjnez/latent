library(lavaan)
library(latent)

samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
full <- hexaco[, selection]
mooc <- full[hexaco$sample == "mooc", ]
dim(mooc)

model.EM <- "FEA =~ hexemfea146 + hexemfea170 + hexemfea74 + hexemfea2
             ANX =~ hexemanx128 + hexemanx8 + hexemanx80 + hexemanx176
             DEP =~ hexemdep62 + hexemdep182 + hexemdep134 + hexemdep158
             SEN =~ hexemsen44 + hexemsen164 + hexemsen20 + hexemsen68"

estimator <- "dwls"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

fit <- lcfa(data = mooc,
            model = model.EM,
            ordered = TRUE,
            estimator = estimator,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            se = "standard",
            control = NULL,
            do.fit = TRUE)

# With lavaan:
fit2 <- lavaan::cfa(data = mooc,
                    model = model.EM,
                    ordered = TRUE,
                    estimator = estimator,
                    meanstructure = meanstructure,
                    std.ov = std.ov,
                    std.lv = std.lv,
                    parameterization = "delta",
                    do.fit = TRUE)

# The two-step polychoric weight estimates differ between the packages.
latent_loss <- latInspect(fit, "loss")["loss", "overall"]
n <- nrow(mooc)
c(latent = latent_loss, lavaan_default = lavaan::lavInspect(fit2, "test")$standard$stat/(n-1))

# latent sums both halves of the correlation matrix and half-weights the
# threshold residuals. Use those weights in lavaan's unique-statistic order.
est <- fit@modelInfo$control_estimator[[1L]]
WLS.V <- diag(c(est$w_thresholds/2, est$W[lower.tri(est$W)]))
fit2_matched <- lavaan::cfa(model.EM, data = mooc, ordered = TRUE,
                            estimator = estimator, meanstructure = meanstructure,
                            std.lv = std.lv, std.ov = std.ov,
                            parameterization = "delta", WLS.V = WLS.V)
lavaan_loss <- lavaan::lavInspect(fit2_matched, "test")$standard$stat/(n-1)
c(latent = latent_loss, lavaan_matched = lavaan_loss)
stopifnot(abs(latent_loss - lavaan_loss) < 1e-6)

diag(lavaan::lavInspect(fit2_matched, "est")$theta)
diag(latInspect(fit, "est", sort = FALSE)$theta)
