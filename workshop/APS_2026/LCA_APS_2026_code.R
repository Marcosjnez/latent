# R code extracted from LCA_APS_2026.qmd
# Updated from the attached LCA_APS_2026_code.R.
# Source of truth: R chunks in the Quarto file, kept in their original order.
# Quarto chunk options such as #| fig-width are preserved as comments.

# ---- Effects of ignoring ----
# CFA example: create scale scores, fit the measurement model, and compare observed and latent correlations.

#| fig-width: 8
#| fig-height: 5

library(lavaan)
library(semPlot)

dat <- HolzingerSwineford1939
dat$vis <- rowMeans(dat[,c("x1","x2","x3")])
dat$tex <- rowMeans(dat[,c("x4","x5","x6")])
dat$spd <- rowMeans(dat[,c("x7","x8","x9")])

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data = HolzingerSwineford1939, std.lv=T)

semPaths(fit)

## average items
round(cor(dat[,c("vis","tex","spd")]),3)

## measurement model
lavInspect(fit, "cor.lv")

# ---- The latent R package ----
# Load the latent package used in the LCA examples.

library(latent)

# ---- LCA example ----
# Categorical-indicator LCA example using the ESS data.

library(collapse)
library(summarytools)
library(stevemisc)
library(ggplot2)
library(patchwork)

# ESS_datLCA <- rio::import("ESS_datLCA.sav")
dat <- sbt(ESS_datLCA, cntry == "FR")
dat <- tidyr::drop_na(dat)

dat[,1:8] <- dapply(dat[,1:8], MARGIN = 2,
                    function(x) carrec(x,"1='Yes';2='No'"))

dat$gndr <- carrec(dat$gndr, "1='Men'; 2 ='Women' ")

freq(dat[,1:2])

freq(dat[,3:4])

freq(dat[,5:6])

freq(dat[,7:8])

fit1 <- lca(data = dat, nclasses = 1:7,
           multinomial = c("contplt", "badge", "sgnptit", "bctprd", 
                           "donprty", "pbldmna", "volunfp","pstplonl") )

fit_ind <- getfit(fit1, digits = 2)
fit_ind

#| fig-width: 8
#| fig-height: 5

plot(fit_ind)

pr2 <- latInspect(fit1[[2]], "profile")
pr2$class
pr2$item[1:3]

pr2$item[4:8]

pr3 <- latInspect(fit1[[3]], "profile")
pr3$class
pr3$item[1:3]

pr3$item[4:8]

#| fig-width: 8
#| fig-height: 5

fit1_3 <- fit1[[3]]

plot(fit1_3)

latInspect(fit1_3, "classes")

# ---- LCA example: Classification Diagnostics ----
# Classification-quality diagnostics for the selected LCA solution.

cl_diag <- lclass_diag(fit1_3, digits = 3, type = "all")

cl_diag$Sum.Posterior

cl_diag$Sum.Mostlikely

cl_diag$AvePP

cl_diag$Entropy

cl_diag

# ---- LCA example: predicted classes ----
# Predicted class assignment and cross-tabulations.

latInspect(fit1_3, "state")

dat$State <- latInspect(fit1_3, "state")
ctable(dat$gndr, dat$State, prop="c")
ctable(dat$gndr, dat$State, prop="r")

# ---- Global model fit evaluation ----
# Likelihood-ratio fit statistics for the selected model.

getfit(fit1_3)[c("L2","dof","pvalue")]

# ---- Residual evaluation in latent ----
# Bivariate residual diagnostics and residual-dependency examples.

res_diag <- lbvr(fit1_3, digits = 3)
res_diag

fit1_3_resid <- lca(data = dat, nclasses = 3,
                    multinomial = c("contplt", "badge", "sgnptit", "bctprd",
                                    "donprty", "pbldmna", "volunfp","pstplonl"),
                    model = list("contplt ~~ volunfp"))

ci_fit1_3_resid <- ci(fit1_3_resid, confidence = 0.95)
ci_fit1_3_resid$table$log_contpltxlog_volunfp

# ---- LCA with continuous indicators ----
# Continuous-indicator LCA example using Big Five scale scores.

bfi <- psych::bfi

vars_rev <- c("A1","C4","C5","E1","E2","O1","O3","O4")

bfi[,vars_rev] <- dapply(bfi[,vars_rev], MARGIN = 2,
                    function(x) carrec(x,"1=6; 2=5; 3=4; 4=3; 5=2; 6=1"))


bfi$Agreeableness <- rowMeans(bfi[,paste0("A",1:5)], na.rm=T)
bfi$Conscientiousness <- rowMeans(bfi[,paste0("C",1:5)], na.rm=T)
bfi$Extraversion <- rowMeans(bfi[,paste0("E",1:5)], na.rm=T)
bfi$Neuroticism <- rowMeans(bfi[,paste0("N",1:5)], na.rm=T)
bfi$Openness <- rowMeans(bfi[,paste0("O",1:5)], na.rm=T)

#head(bfi)

summarytools::descr(bfi[,29:33])

#| fig-width: 8
#| fig-height: 5

ggplot(bfi, aes(x=Agreeableness))+geom_density()+
  ggplot(bfi, aes(x=Conscientiousness))+geom_density()+
  ggplot(bfi, aes(x=Extraversion))+geom_density()+
  ggplot(bfi, aes(x=Neuroticism))+geom_density()+
  ggplot(bfi, aes(x=Openness))+geom_density()

#| output: false

set.seed(1987)

fit2 <- lca(data = bfi, nclasses = 1:7,
           gaussian = c("Agreeableness", "Conscientiousness", 
                        "Extraversion", "Neuroticism", "Openness" ) )

fit_ind_cont <- getfit(fit2, digits = 2)
fit_ind_cont

#| fig-width: 8
#| fig-height: 5

plot(fit_ind_cont)

latInspect(fit2[[3]], "classes")
latInspect(fit2[[4]], "classes")
latInspect(fit2[[5]], "classes")

fit2_4 <- fit2[[4]]

pr4 <- latInspect(fit2_4, "profile")
pr4$class
pr4$item[1:2]

pr4$item[3:5]

#| fig-width: 8
#| fig-height: 5

plot(fit2_4)

# ---- LCA example ----
# Categorical-indicator LCA example using the ESS data.

latInspect(fit2_4, "classes")

# ---- LCA example: Classification Diagnostics ----
# Classification-quality diagnostics for the selected LCA solution.

cl_diag_cont <- lclass_diag(fit2_4, type = "all", digits = 3)

cl_diag_cont$Sum.Posterior

cl_diag_cont$Sum.Mostlikely

cl_diag_cont$AvePP

cl_diag_cont$Entropy

cl_diag_cont

# ---- LCA example: predicted classes ----
# Predicted class assignment and cross-tabulations.

latInspect(fit2_4, "state")

bfi$Gender <- carrec(bfi$gender, "1='Man';2='Woman' ")
bfi$State <- as.factor(latInspect(fit2_4, "state"))
ctable(bfi$Gender, bfi$State, prop="c")
ctable(bfi$Gender, bfi$State, prop="r")

# ---- LCA with continuous indicators ----
# Continuous-indicator LCA example using Big Five scale scores.

#| fig-width: 8
#| fig-height: 5

ggplot(bfi, aes(x=Agreeableness, group = State))+
  geom_density(aes(color = State))+
  ggplot(bfi, aes(x=Conscientiousness, group = State))+
  geom_density(aes(color = State))+
  ggplot(bfi, aes(x=Extraversion, group = State))+
  geom_density(aes(color = State))+
  ggplot(bfi, aes(x=Neuroticism, group = State))+
  geom_density(aes(color = State))+
  ggplot(bfi, aes(x=Openness, group = State))+
  geom_density(aes(color = State))

# ---- Residual evaluation in latent ----
# Bivariate residual diagnostics and residual-dependency examples.

res_diag2 <- lbvr(fit2_4, digits = 3)
res_diag2

fit2_4_resid <- lca(data = bfi, nclasses = 4,
            gaussian = c("Agreeableness", "Conscientiousness",
                         "Extraversion", "Neuroticism", "Openness" ),
            model = list("Agreeableness ~~ Extraversion"))

getfit(fit2_4)
getfit(fit2_4_resid)

# ---- LCA with continuous and categorical items ----
# Mixed LCA example combining continuous and categorical indicators.

# dat3 <- rio::import("hs_academic.sav")
dat3 <- hs_academic

vars_rev <- c("hm","he","voc","nocol")

dat3[,vars_rev] <- dapply(dat3[,vars_rev], MARGIN = 2,
                    function(x) carrec(x,"0 = 'No'; 1 = 'Yes' "))

head(dat3)

freq(dat3[,1:4])

summarytools::descr(dat3[,5:8])

fit3 <- lca(data = dat3, nclasses = 1:7,
           gaussian = c("ach9","ach10","ach11","ach12"),
           multinomial = c("hm","he","voc","nocol"))

fit_ind3 <- getfit(fit3, digits = 2)
fit_ind3

plot(fit_ind3)

fit3_2 <- fit3[[2]]

pr2 <- latInspect(fit3_2, "profile")
pr2$class
pr2$item[1:4]

pr2$item[5:8]

plots_f3 <- plot(fit3_2)

plots_f3[[1]]/plots_f3[[2]]

cl_diag_mix <- lclass_diag(fit3_2, type = "all", digits = 3)

cl_diag_mix$Sum.Posterior

cl_diag_mix$Sum.Mostlikely

cl_diag_mix$AvePP

cl_diag_mix$Entropy

cl_diag_cont

# ---- Residual evaluation in latent ----
# Bivariate residual diagnostics and residual-dependency examples.

res_diag3 <- lbvr(fit3_2, digits = 3)
res_diag3

# ---- Adding predictors - latent ----
# Add predictors to the latent class model and inspect their effects.

fit4 <- lca(data = dat, nclasses = 3,
           multinomial = c("contplt", "badge", "sgnptit", "bctprd", 
                           "donprty", "pbldmna", "volunfp","pstplonl") )

fit4_step2 <- lca(data = dat,
                  nclasses = 3,
                  multinomial = c("contplt", "badge", 
                                  "sgnptit", "bctprd", 
                                  "donprty", "pbldmna", 
                                  "volunfp","pstplonl"),
                 X = c("gndr", "agea"),
                 model = fit4)

latInspect(fit4_step2, "coefs")

se_step2 <- se(fit4_step2)
se_step2$table$beta

# Compute 95% confidence intervals for predictor coefficients.
ci_step2 <- ci(fit4_step2, confidence = 0.95)
ci_step2$table$beta

# Note: this chunk was marked eval=FALSE in the .qmd; run manually if needed.
x <- plot_coeffs(fit4_step2,
                 what = "OR",        # Odds ratios
                 effects = "coding", # Effects-coding parametrization
                 confidence = 0.95,  # Confidence level
                 xlim = c(0, 3),     # Plot limits for x-axis
                 mfrow = c(2, 2))    # Grid of viewer panel

head(predict(fit4_step2))

predict(fit4_step2, 
        new = rbind(c(0,50),
                    c(1,50)))

predict(fit4_step2, 
        new = rbind(c(0,45),
                    c(0,50),
                    c(1,45),
                    c(1,50)))

# ---- Assessing convergence ----
# Examples for checking convergence, log-likelihood, and gradients.

set.seed(2026)
fit <- lca(data = empathy,
           nclasses = 2L,
           gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
           control = list(maxit = 20))

latInspect(fit, what = "loglik")
latInspect(fit, what = "convergence")
latInspect(fit, what = "gradient")

set.seed(2026)
fit <- lca(data = cancer,
           nclasses = 3L,
           gaussian = c("Age", "WeightIndex", "SystolicBloodPressure",
                        "DiastolicBloodPressure"),
           multinomial = c("PerformanceRating", "CardiovascularDiseaseHistory"),
           control = list(rstarts = 1L))

latInspect(fit, what = "loglik")
latInspect(fit, what = "convergence")
latInspect(fit, what = "gradient")

control = list(opt     = "lbfgs", # newton
               maxit   = 1000L,   # Maximum number of iterations
               rstarts = 16L,     # Number of random starts
               cores   = 1L,      # Number or parallel runs
               eps     = 1e-05)   # Convergence criteria (gradient norm)

# ---- Penalties ----
# Examples of penalized estimation and penalty defaults.

set.seed(59) # Pick a seed that creates the problem
fit_empathy <- lca(data = empathy,
           nclasses = 4L,
           gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
           penalties = FALSE) # Deactivate penalized-ML, just ML

latInspect(fit_empathy, what = "loglik")
latInspect(fit_empathy, what = "convergence")
latInspect(fit_empathy, what = "profile") # What's wrong in this profile?

set.seed(59) # Pick a seed that creates the problem
fit_empathy <- lca(data = empathy,
           nclasses = 4L,
           gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
           penalties = TRUE)

latInspect(fit_empathy, what = "loglik")
latInspect(fit_empathy, what = "convergence")
latInspect(fit_empathy, what = "profile")

penalty_defaults <- list(
  class = list(alpha = 1), # Penalty for class probabilities
  prob  = list(alpha = 1), # Penalty for item probabilities (multinomial items)
  var   = list(alpha = 1), # Penalty for item variances (gaussian items)
  Sigma = list(alpha = 1)  # Penalty for covariance matrices (gaussian items)
)

# Note: this chunk was marked eval=FALSE in the .qmd; run manually if needed.
latInspect(fit_empathy, what = "profile") # Before
fit_empathy <- lca(data = empathy,
                   nclasses = 4L,
                   gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
                   penalties = list(var = list(alpha = 100))) # Extreme penalty
latInspect(fit_empathy, what = "profile") # After
