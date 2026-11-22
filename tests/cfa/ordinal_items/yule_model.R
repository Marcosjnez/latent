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

model.EM <- "FEA =~ hexemfea146 + hexemfea170 + hexemfea74 + hexemfea2
             ANX =~ hexemanx128 + hexemanx8 + hexemanx80 + hexemanx176
             DEP =~ hexemdep62 + hexemdep182 + hexemdep134 + hexemdep158
             SEN =~ hexemsen44 + hexemsen164 + hexemsen20 + hexemsen68"

estimator <- "dwls"
std.ov <- FALSE
std.lv <- FALSE
meanstructure <- TRUE
likelihood <- "normal"

fit <- lcfa(model = model.EM,
            data = mooc,
            ordered = "yule",
            estimator = estimator,
            positive = FALSE,
            penalties = FALSE,
            std.ov = std.ov,
            std.lv = std.lv,
            meanstructure = meanstructure,
            se = "standard",
            do.fit = TRUE,
            control = NULL)
latInspect(fit, "loss")
fit@Optim$iterations
fit@Optim$convergence
fit@timing

latInspect(fit, what = "lambda")
# This Yule example uses an approximate inverse for the first-stage Hessian
# when propagating the sampling uncertainty into the factor parameters.
fit@Optim$SE$se
