# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 11/06/2026

#### gss82 ####

library(latent)
gss82$EDUCR <- as.integer(gss82$EDUCR)-1L
indicators <- c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT")

set.seed(2026)
fit1 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = indicators,
            penalties = list(class = list(alpha=1),
                             prob  = list(alpha=0)),
            do.fit = TRUE)

latInspect(fit1, what = "loglik")
# loglik: -3891.472 # penalized_loglik: -3896.468
latInspect(fit1, what = "convergence")
latInspect(fit1, what = "profile", digits = 3)

set.seed(2026)
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
# loglik: -3741.225 # penalized_loglik: -3746.363
latInspect(fit2, what = "convergence")

# check that the measurement model was fixed:
all.equal(fit1@parameters[-1], fit2@parameters[-1])

# Inspect model objects:
latInspect(fit2, what = "coefs", digits = 5)
latInspect(fit2, what = "profile", digits = 3)

# Standard errors:
SE2 <- se(fit2, type = "standard", digits = 4)
SE2$se

#### empathy ####

library(latent)

set.seed(2026)
fit1 <- lca(data = empathy,
            nclasses = 4L,
            gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
            do.fit = TRUE)

latInspect(fit1, what = "loglik")
# loglik: -1841.336 # penalized_loglik: -1844.333

latInspect(fit1, what = "convergence")

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
latInspect(fit2, what = "coefs", digits = 5)
latInspect(fit2, what = "profile", digits = 3)

# Standard errors:
SE2 <- se(fit2, type = "standard", digits = 4)
SE2$se
