# Author: Marcos Jimenez
# email: m.j.jimenezhenriquez@vu.nl
# Modification date: 15/06/2026

library(latent)

#### Step 1: Measurement model ####

gss82$EDUCR <- as.integer(gss82$EDUCR) - 1L
indicators <- c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT")

set.seed(2027)
fit1 <- lca(data = gss82,
            nclasses = 3L,
            multinomial = indicators,
            penalties = list(class = list(alpha = 1),
                             prob  = list(alpha = 0)),
            do.fit = TRUE)

post <- latInspect(fit1, what = "posterior")
K <- ncol(post)
N <- nrow(post)

class_error_prop <- latInspect(fit1, what = "classification")$class_error_prop

#### Step 2: BCH weights with proportional classification ####

W_bch_prop <- bch_weights(post,
                          class_error = class_error_prop,
                          type = "proportional")

# Useful diagnostics
rowSums(W_bch_prop)
summary(as.vector(W_bch_prop))
any(W_bch_prop < 0)

#### Step 3: Expanded BCH data ####

gss82_bch_prop <- gss82[rep(seq_len(N), each = K), , drop = FALSE]
row.names(gss82_bch_prop) <- NULL

gss82_bch_prop$states <- factor(rep(seq_len(K), times = N),
                                levels = seq_len(K))

weights_bch_prop <- as.vector(t(W_bch_prop))

eps <- 1e-12

perfect <- matrix(eps, nrow = K, ncol = K)
diag(perfect) <- 1 - eps * (K - 1L)

log_perfect <- log(perfect)
log_perfect <- apply(log_perfect, MARGIN = 1L,
                     FUN = \(x) x - x[1L])

set.seed(2027)

fit_bch_prop <- lca(data = gss82_bch_prop,
                    nclasses = K,
                    multinomial = "states",
                    X = c("RACE", "SEX", "EDUCR", "AGE"),
                    model = list(log_states = log_perfect),
                    penalties = list(class = list(alpha = 1),
                                     prob  = list(alpha = 0)),
                    weights = weights_bch_prop,
                    do.fit = TRUE)

latInspect(fit_bch_prop, what = "loglik")
latInspect(fit_bch_prop, what = "convergence")
latInspect(fit_bch_prop, what = "profile")
latInspect(fit_bch_prop, what = "coefs")
