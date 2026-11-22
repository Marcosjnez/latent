library(latent)
fit <- lcfa(HolzingerSwineford1939[, paste("x", 1:9, sep = "")],
            model = 'visual =~ x1 + x2 + x3
                     textual =~ x4 + x5 + x6
                     speed =~ x7 + x8 + x9')
fs <- latscores(fit, method = "Green")
fs@parameters$scores
fs@Optim$SE$table_se$scores
latInspect(fs, "scores")               # Display order and signs.
se(fs)$table_se                        # Native order and signs.
vcov(fs, person = 1)                   # Small error covariance for person 1.
fs_new <- latscores(fit,
                    data = HolzingerSwineford1939[, paste("x", 1:9, sep = "")],
                    method = "regression")
