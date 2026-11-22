library(latent)
library(lavaan)

samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
full <- hexaco[, selection]

mooc <- full[hexaco$sample == "mooc", ]
dim(mooc)
mooc$school <- rep(c("s1", "s2"), times = c(2000L, nrow(mooc) - 2000L))

model.EM <- "FEA =~ hexemfea146 + hexemfea170 + hexemfea74 + hexemfea2
             ANX =~ hexemanx128 + hexemanx8 + hexemanx80 + hexemanx176
             DEP =~ hexemdep62 + hexemdep182 + hexemdep134 + hexemdep158
             SEN =~ hexemsen44 + hexemsen164 + hexemsen20 + hexemsen68"

fit <- lcfa(model = model.EM,
            data = mooc,
            ordered = TRUE,
            estimator = "dwls",
            group = "school",
            se = "standard",
            do.fit = TRUE,
            control = NULL)
loss <- latInspect(fit, "loss")["loss", ]

# With lavaan, use the same marker indicators and delta parameterization.
fit2 <- lavaan::cfa(model = model.EM, data = mooc, ordered = TRUE,
                    estimator = "dwls", group = "school", std.lv = FALSE,
                    meanstructure = TRUE, parameterization = "delta")
n <- as.numeric(lavaan::lavInspect(fit2, "nobs"))
default_loss <- lavaan::lavInspect(fit2, "test")$standard$stat.group/(n-1)
c(latent = unname(loss["overall"]), lavaan_default = weighted.mean(default_loss, n))

WLS.V <- lapply(fit@modelInfo$control_estimator, function(est) {
  diag(c(est$w_thresholds/2, est$W[lower.tri(est$W)]))
})
fit2_matched <- lavaan::cfa(model.EM, data = mooc, ordered = TRUE,
                            estimator = "dwls", group = "school", std.lv = FALSE,
                            meanstructure = TRUE, parameterization = "delta",
                            WLS.V = WLS.V)
group_loss <- lavaan::lavInspect(fit2_matched, "test")$standard$stat.group/(n-1)
# latent reports each group's contribution weighted by its sample fraction.
comparison <- cbind(latent = loss[-length(loss)], lavaan_matched = group_loss*n/sum(n))
comparison
stopifnot(max(abs(comparison[, 1L] - comparison[, 2L])) < 1e-6)

# Inspect sampling uncertainty using the group-specific parameter structures.
blocks <- with(fit@dataList$data_param, c(lambda_group, theta_group, psi_group))
latent_se <- se(fit, parameters = fit@modelInfo$trans[blocks], digits = NULL)$table_se
lavaan_se <- lavaan::lavInspect(fit2, "se")
latent_se[c("lambda.s1", "theta.s1", "psi.s1")]
lavaan_se$s1[c("lambda", "theta", "psi")]

# Compare correlation weights separately from thresholds, in the same pair order.
W <- lavaan::lavInspect(fit2, "wls.v")
lapply(seq_along(W), function(g) {
  native <- fit@modelInfo$control_estimator[[g]]$W
  weights <- diag(W[[g]])
  pairs <- grepl("~~", names(weights), fixed = TRUE)
  cbind(latent = native[lower.tri(native)], lavaan = weights[pairs])
})
