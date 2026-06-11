# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 11/06/2026

# This code runs the step 3 method of LG with options:
# Analysis: covariates
# Classification: Modal
# Adjustment: ML

library(latent)

#### Measurement model ####

gss82$EDUCR <- as.integer(gss82$EDUCR)-1L
indicators <- c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT")

set.seed(2026)
fit1 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = indicators,
            do.fit = TRUE)

latInspect(fit1, what = "loglik")
# loglik: -3891.252 # penalized_loglik: -3892.478
latInspect(fit1, what = "convergence")
latInspect(fit1, what = "profile")

#### Modal assignment ####

# post <- latInspect(fit1, what = "posterior")
gss82$states <- latInspect(fit1, what = "state")

fixed <- lclass_diag(fit1)$Mostlikely.Class # Classification table
logfixed <- log(fixed)
logfixed <- apply(logfixed, MARGIN = 1L, FUN = \(x) x-x[1])
# apply(logfixed, 2, soft, a=1)

#### Fitting the covariates model using the states ####

# RACE and SEX are treated as nominal and EDUCR and AGE as continuous
fit2 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = "states",
            X = c("RACE", "SEX", "EDUCR", "AGE"),
            model = list(log_states = logfixed),
            do.fit = TRUE)

latInspect(fit2, what = "loglik")
# loglik: -1581.264 # penalized_loglik: -1582.497
latInspect(fit2, what = "convergence")

latInspect(fit2, what = "profile") # CLUSTER SIZE DIFFERS FROM LG, WHY?
latInspect(fit2, what = "coefs")

# NO VALID STANDARD ERRORS BECAUSE VAR(logfixed) IS NOT CALCULATED IN fit1
# Solution:
# (1) Calculate logfixed as transformed parameters in fit1
# (2) Provide the full fit1 with logfixed in model
# (3) Automatically propagate SE estimation for transformed parameters
