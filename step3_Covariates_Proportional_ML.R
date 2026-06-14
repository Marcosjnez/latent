# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 14/06/2026

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

N <- nrow(gss82)
gss82 <- gss82[rep(seq_len(N), each = 3L), ] # Expand the dataset
gss82$states <- factor(rep(seq_len(3), times = N), levels = seq_len(3))
posterior <- latInspect(fit1, what = "posterior")
weights <- as.vector(t(posterior))

# Proportional classification-error matrix:
# rows = true latent class C
# columns = assigned/pseudo class W
fixed_prop <- crossprod(posterior, posterior)
fixed_prop <- sweep(fixed_prop, 1L, colSums(posterior), "/")
# Convert to the log-parameterization used by latent
logfixed_prop <- log(fixed_prop)
logfixed_prop <- apply(logfixed_prop, MARGIN = 1L, FUN = \(x) x - x[1L])

#### Step 3: Fitting the covariate model using the states ####

# RACE and SEX are treated as nominal and EDUCR and AGE as continuous
set.seed(2027)
fit2 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = "states",
            X = c("RACE", "SEX", "EDUCR", "AGE"),
            model = list(log_states = logfixed_prop),
            penalties = list(class = list(alpha=1),
                             prob  = list(alpha=0)),
            weights = weights,
            do.fit = TRUE)

latInspect(fit2, what = "loglik")
# loglik: -1479.613 # penalized_loglik: -1482.190
latInspect(fit2, what = "convergence")

latInspect(fit2, what = "profile")
latInspect(fit2, what = "coefs")

# Standard errors:
SE <- se(fit2, type = "standard")
SE$se
