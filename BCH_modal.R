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

latInspect(fit1, what = "loglik")
# loglik: -3891.252 # penalized_loglik: -3892.478
post <- latInspect(fit1, what = "posterior")
K <- ncol(post)
N <- nrow(post)

class_error_modal <- latInspect(fit1, what = "classification")$class_error_modal

#### Step 2: BCH weights with modal classification ####

W_bch_modal <- bch_weights(post,
                           class_error = class_error_modal,
                           type = "modal")

# Useful diagnostics
rowSums(W_bch_modal)
summary(as.vector(W_bch_modal))
any(W_bch_modal < 0)

#### Step 3: Expanded BCH data ####

gss82_bch_modal <- gss82[rep(seq_len(N), each = K), , drop = FALSE]
row.names(gss82_bch_modal) <- NULL

gss82_bch_modal$states <- factor(rep(seq_len(K), times = N),
                                 levels = seq_len(K))

weights_bch_modal <- as.vector(t(W_bch_modal))

#### Perfect measurement matrix for states ####

eps <- 1

perfect <- matrix(eps, nrow = K, ncol = K)
diag(perfect) <- 1 - eps * (K - 1L)

log_perfect <- log(perfect)
log_perfect <- apply(log_perfect, MARGIN = 1L,
                     FUN = \(x) x - x[1L])

class_error_modal <- latInspect(fit1, what = "classification")$class_error_modal
# Convert to the log-parameterization used by latent
log_class_error_modal <- log(class_error_modal)
log_class_error_modal <- apply(log_class_error_modal, MARGIN = 1L,
                               FUN = \(x) x - x[1L])
t(apply(log_class_error_modal, 2, soft, a=1)) / class_error_modal

set.seed(2027)

fit_bch_modal <- lca(data = gss82_bch_modal,
                     nclasses = K,
                     multinomial = "states",
                     X = c("RACE", "SEX", "EDUCR", "AGE"),
                     model = list(log_states = log_class_error_modal),
                     penalties = list(class = list(alpha = 1),
                                      prob  = list(alpha = 0)),
                     weights = weights_bch_modal,
                     do.fit = TRUE)

latInspect(fit_bch_modal, what = "loglik")
latInspect(fit_bch_modal, what = "convergence")
latInspect(fit_bch_modal, what = "profile")
latInspect(fit_bch_modal, what = "coefs")
