# GET THE SAME LOGLIK THAN LATENTGOLD FOR STRUCTURAL FIT

library(latent)
set.seed(2026)

gss82$EDUCR <- as.integer(gss82$EDUCR)-1L
fit <- lca(data = gss82,
           nclasses = 3L,
           multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
           covariates = c("RACE", "SEX", "EDUCR", "AGE"),
           penalties = list(class = list(alpha=1),
                            prob  = list(alpha=0)),
           do.fit = TRUE)
latInspect(fit, what = "loglik")

gss82 <- na.omit(gss82[, c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT",
                           "RACE", "SEX", "EDUCR", "AGE")])
measurement <- lca(data = gss82,
                   nclasses = 3L,
                   multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
                   penalties = list(class = list(alpha=1),
                                    prob  = list(alpha=0)),
                   model = fit$measurement,
                   control = list(free_beta = FALSE),
                   do.fit = TRUE)
latInspect(measurement, what = "loglik")
latInspect(measurement, what = "convergence")

latInspect(fit$structural, what = "loglik") -
  latInspect(measurement, what = "loglik")
