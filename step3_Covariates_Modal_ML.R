# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/06/2026

# This code runs the step 3 method of LG with options:
# Analysis: Covariates
# Classification: Modal
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

#### Step 2: Modal assignment ####

gss82$states <- latInspect(fit1, what = "state")
# Modal classification-error matrix:
# rows = true latent class C
# columns = assigned/pseudo class W
class_error_modal <- latInspect(fit1, what = "classification")$class_error_modal
# Convert to the log-parameterization used by latent
log_class_error_modal <- log(class_error_modal)
log_class_error_modal <- apply(log_class_error_modal, MARGIN = 1L,
                               FUN = \(x) x - x[1L])
# t(apply(log_class_error_modal, 2, soft, a=1)) / class_error_modal

#### Step 3: Fitting the covariate model using the states ####

# RACE and SEX are treated as nominal and EDUCR and AGE as continuous
fit2 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = "states",
            X = c("RACE", "SEX", "EDUCR", "AGE"),
            model = list(log_states = log_class_error_modal),
            penalties = list(class = list(alpha=1),
                             prob  = list(alpha=0)),
            adjustment = "none",
            do.fit = TRUE)

latInspect(fit2, what = "loglik")
# loglik: -1483.887 # penalized_loglik: -1485.214
latInspect(fit2, what = "convergence")

latInspect(fit2, what = "profile")
latInspect(fit2, what = "coefs")

# NO VALID STANDARD ERRORS BECAUSE VAR(log_class_error_modal) IS NOT CALCULATED IN fit1
# Solution:
# (1) Calculate log_class_error_modal as transformed parameters in fit1
# (2) Provide the full fit1 with log_class_error_modal in model
# (3) Automatically propagate SE estimation for transformed parameters
