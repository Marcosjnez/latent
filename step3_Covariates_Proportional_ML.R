# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/06/2026

# This code runs the step 3 method of LG with options:
# Analysis: Covariates
# Classification: Proportional
# Adjustment: ML

library(latent)

#### Step 1: Measurement model ####

gss82$EDUCR <- as.integer(gss82$EDUCR)-1L
indicators <- c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT")

set.seed(2027)
fit1 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = indicators,
            penalties = list(class = list(alpha=1),
                             prob  = list(alpha=0)),
            do.fit = TRUE)

latInspect(fit1, what = "loglik")
# loglik: -3891.252 # penalized_loglik: -3892.478
latInspect(fit1, what = "convergence")
latInspect(fit1, what = "profile")

#### Step 2: Compute the weights ####

gss82$states <- latInspect(fit1, what = "state")
weights <- latInspect(fit1, what = "posterior")

#### Step 3: Fitting the covariate model using the states ####

# RACE and SEX are treated as nominal and EDUCR and AGE as continuous
set.seed(2027)
weights <- latInspect(fit1, what = "posterior")
fit2 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = "states",
            X = c("RACE", "SEX", "EDUCR", "AGE"),
            penalties = list(class = list(alpha=1),
                             prob  = list(alpha=1)),
            weights = weights,
            do.fit = TRUE)

latInspect(fit2, what = "loglik")
# loglik: -1479.806 # penalized_loglik: -1483.852
latInspect(fit2, what = "convergence")

latInspect(fit2, what = "profile")
latInspect(fit2, what = "coefs")

# Standard errors:
SE <- se(fit2, type = "standard")
SE$se
