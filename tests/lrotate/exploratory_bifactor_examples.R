library(latent)
library(GPArotation)

set.seed(2026)
ulambda <- rorth(9, 4)

# Bi-geomin:

rot <- bigeominQ(ulambda, delta = 0.01, randomStarts = 10L)
fit <- lrotate(lambda = ulambda, projection = "oblq",
               rotation = list(geomin = list(factors = 2:4, epsilon = 0.01)),
               control = list(rstarts = 10L))
fit@Optim$iterations
fit@Optim$f
rotation_differences(fit, rot$loadings, rot$Phi)

# Bi-quartimin:

rot <- bifactorQ(ulambda, randomStarts = 10L)
fit <- lrotate(lambda = ulambda, projection = "oblq",
               rotation = list(oblimin = list(factors = 2:4)),
               control = list(rstarts = 10L))
fit@Optim$iterations
fit@Optim$f
rotation_differences(fit, rot$loadings, rot$Phi)
