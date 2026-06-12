# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 12/06/2026

# This code runs the step 3 method of LG with options:
# Analysis: Covariates
# Classification: Proportional
# Adjustment: ML

library(latent)

#### Measurement model ####

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

#### Compute the weights ####

weights <- latInspect(fit1, what = "posterior")

#### Fitting the covariates model using the states ####

# RACE and SEX are treated as nominal and EDUCR and AGE as continuous
set.seed(2027)
fit2 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = indicators,
            X = c("RACE", "SEX", "EDUCR", "AGE"),
            model = fit1,
            penalties = list(class = list(alpha=1),
                             prob  = list(alpha=0)),
            weights = weights,
            do.fit = TRUE)

latInspect(fit2, what = "loglik")
# loglik: -3739.894 # penalized_loglik: -3741.265
latInspect(fit2, what = "convergence")

latInspect(fit2, what = "profile") # CLUSTER SIZE DIFFERS FROM LG, WHY?
latInspect(fit2, what = "coefs")

# I THINK STANDARD ERRORS SHOULD BE VALID

# Standard errors:
SE <- se(fit2, type = "standard", digits = 4)
SE$se
