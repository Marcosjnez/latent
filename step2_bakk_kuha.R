# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 13/06/2026

library(latent)

#### gss82 ####

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
# loglik: -3891.472 # penalized_loglik: -3896.468
latInspect(fit1, what = "convergence")
latInspect(fit1, what = "profile")

#### Step 2: Fitting the covariate model fixing the measurement part ####

fit2 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = indicators,
            X = c("RACE", "SEX", "EDUCR", "AGE"),
            model = fit1,
            # model = list("UNDERSTA ~~ COOPERAT
            #               UNDERSTA ~~ ACCURACY"),
            penalties = list(class = list(alpha=1),
                             prob  = list(alpha=0)),
            do.fit = TRUE)

latInspect(fit2, what = "loglik")
# loglik: -3739.894 # penalized_loglik: -3741.265
latInspect(fit2, what = "convergence")

# check that the measurement model was fixed:
all.equal(fit1@parameters[-1], fit2@parameters[-1])

# Inspect model objects:
latInspect(fit2, what = "coefs")
latInspect(fit2, what = "profile")

# Standard errors:
SE <- se(fit2, type = "standard")
SE$se

#### empathy ####

#### Step 1: Measurement model ####

set.seed(2026)
fit1 <- lca(data = empathy,
            nclasses = 4L,
            gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
            do.fit = TRUE)

latInspect(fit1, what = "loglik")
# loglik: -1841.336 # penalized_loglik: -1844.333

latInspect(fit1, what = "convergence")

#### Step 2: Fitting the covariate model fixing the measurement part ####

set.seed(2026)
fit2 <- lca(data = empathy,
            nclasses = 4L,
            gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
            X = c("pt1", "pt2", "pt3", "pt4"),
            model = fit1,
            do.fit = TRUE)

latInspect(fit2, what = "loglik")
# loglik: -1747.135 # penalized_loglik: -1750.566

latInspect(fit2, what = "convergence")

# check that the measurement model was fixed:
all.equal(fit1@parameters[-1], fit2@parameters[-1])

# Inspect model objects:
latInspect(fit2, what = "coefs")
latInspect(fit2, what = "profile")

# Standard errors:
SE <- se(fit2, type = "standard")
SE$se
