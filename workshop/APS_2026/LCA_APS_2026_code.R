# R code extracted from LCA_APS_2026.qmd.

# ---- CFA example ----
# Load CFA packages, create scale scores, fit the CFA model, and plot the path diagram.
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

# Compare correlations among averaged item scales.
## average items
round(cor(dat[,c("vis","tex","spd")]),3)

# Inspect correlations among the latent factors.
## measurement model
lavInspect(fit, "cor.lv")


# ---- latent package ----
# Load the latent package for latent class analysis.
library(latent)


# ---- Categorical LCA example ----
# Load helper packages, import ESS data, keep France, recode variables, and inspect frequencies.
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

# Inspect frequencies for the next pair of ESS indicators.
freq(dat[,3:4])

# Inspect frequencies for the next pair of ESS indicators.
freq(dat[,5:6])

# Inspect frequencies for the final pair of ESS indicators.
freq(dat[,7:8])

# Fit LCA models with 1 to 7 classes using categorical indicators.
fit1 <- lca(data = dat, nclasses = 1:7,
           multinomial = c("contplt", "badge", "sgnptit", "bctprd",
                           "donprty", "pbldmna", "volunfp","pstplonl") )

# Extract and display fit indices for the categorical LCA models.
fit_ind <- getfit(fit1, digits = 2)
fit_ind

# Plot model-fit indices across class solutions.
plot(fit_ind)

# Inspect the 2-class profile and the first three items.
pr2 <- latInspect(fit1[[2]], "profile")
pr2$class
pr2$item[1:3]

# Inspect remaining item probabilities for the 2-class model.
pr2$item[4:8]

# Inspect the 3-class profile and the first three items.
pr3 <- latInspect(fit1[[3]], "profile")
pr3$class
pr3$item[1:3]

# Inspect remaining item probabilities for the 3-class model.
pr3$item[4:8]

# Select the 3-class model and plot its class profiles.
fit1_3 <- fit1[[3]]

plot(fit1_3)

# Inspect class proportions for the selected 3-class model.
latInspect(fit1_3, "classes")


# ---- Categorical LCA classification diagnostics ----
# Compute classification diagnostics for the selected categorical LCA model.
cl_diag <- lclass_diag(fit1_3, digits = 3, type = "all")

# Show posterior class totals.
cl_diag$Sum.Posterior

# Show most-likely class totals.
cl_diag$Sum.Mostlikely

# Show average posterior probabilities.
cl_diag$AvePP

# Show entropy-based classification quality.
cl_diag$Entropy

# Display all classification diagnostics.
cl_diag


# ---- Categorical LCA predicted classes and fit ----
# Inspect predicted class membership.
latInspect(fit1_3, "state")

# Add predicted classes to the data and cross-tabulate them by gender.
dat$State <- latInspect(fit1_3, "state")
ctable(dat$gndr, dat$State, prop="c")
ctable(dat$gndr, dat$State, prop="r")

# Check likelihood-ratio fit statistics for the selected model.
getfit(fit1_3)[c("L2","dof","pvalue")]

# Compute bivariate residual diagnostics for the selected model.
res_diag <- lbvr(fit1_3, digits = 3)
res_diag

# Refit the LCA model with 3 classes and the residual dependency between
# `contplt` and `volunfp`.
fit1_3_resid <- lca(data = dat, nclasses = 3,
                    multinomial = c("contplt", "badge", "sgnptit", "bctprd",
                                    "donprty", "pbldmna", "volunfp","pstplonl"),
                    model = list("contplt ~~ volunfp"))

# ---- LCA with continuous indicators ----
# Prepare BFI data, reverse-code items, and compute Big Five scale scores.
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

# Describe the Big Five scale scores.
summarytools::descr(bfi[,29:33])

# Plot density curves for the Big Five scale scores.
ggplot(bfi, aes(x=Agreeableness))+geom_density()+
  ggplot(bfi, aes(x=Conscientiousness))+geom_density()+
  ggplot(bfi, aes(x=Extraversion))+geom_density()+
  ggplot(bfi, aes(x=Neuroticism))+geom_density()+
  ggplot(bfi, aes(x=Openness))+geom_density()

# Fit continuous-indicator LCA models with 1 to 7 classes.
set.seed(1987)

fit2 <- lca(data = bfi, nclasses = 1:7,
           gaussian = c("Agreeableness", "Conscientiousness",
                        "Extraversion", "Neuroticism", "Openness" ) )

# Extract and display fit indices for the continuous LCA models.
fit_ind_cont <- getfit(fit2, digits = 2)
fit_ind_cont

# Plot model-fit indices across continuous class solutions.
plot(fit_ind_cont)

# Compare class proportions for candidate continuous models.
latInspect(fit2[[3]], "classes")
latInspect(fit2[[4]], "classes")
latInspect(fit2[[5]], "classes")

# Select the 4-class continuous model and inspect profile output.
fit2_4 <- fit2[[4]]

pr4 <- latInspect(fit2_4, "profile")
pr4$class
pr4$item[1:2]

# Inspect remaining profile output for the 4-class model.
pr4$item[3:5]

# Plot profiles for the selected continuous model.
plot(fit2_4)

# Inspect class proportions for the selected continuous model.
latInspect(fit2_4, "classes")


# ---- Continuous LCA classification diagnostics ----
# Compute classification diagnostics for the continuous LCA model.
cl_diag_cont <- lclass_diag(fit2_4, type = "all", digits = 3)

# Show posterior class totals.
cl_diag_cont$Sum.Posterior

# Show most-likely class totals.
cl_diag_cont$Sum.Mostlikely

# Show average posterior probabilities.
cl_diag_cont$AvePP

# Show entropy-based classification quality.
cl_diag_cont$Entropy

# Display all classification diagnostics.
cl_diag_cont


# ---- Continuous LCA predicted classes and residuals ----
# Inspect predicted class membership.
latInspect(fit2_4, "state")

# Add predicted classes to the BFI data and cross-tabulate them by gender.
bfi$Gender <- carrec(bfi$gender, "1='Man';2='Woman' ")
bfi$State <- as.factor(latInspect(fit2_4, "state"))
ctable(bfi$Gender, bfi$State, prop="c")
ctable(bfi$Gender, bfi$State, prop="r")

# Plot score distributions by predicted class.
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

# Compute bivariate residual diagnostics for the continuous model.
res_diag2 <- lbvr(fit2_4, digits = 3)
res_diag2

# Refit the model modeling the residual covariance between Agreeableness and
# Extraversion.
fit2_4_resid <- lca(data = bfi, nclasses = 4,
            gaussian = c("Agreeableness", "Conscientiousness",
                         "Extraversion", "Neuroticism", "Openness" ),
            model = list("Agreeableness ~~ Extraversion"))

getfit(fit2_4)
getfit(fit2_4_resid)

# ci_fit2_4_resid <- ci(fit2_4_resid, confidence = 0.95)
# ci_fit2_4_resid$table

# ---- LCA with continuous and categorical indicators ----
# Import mixed-indicator data and recode categorical items.
# dat3 <- rio::import("hs_academic.sav")
dat3 <- hs_academic

vars_rev <- c("hm","he","voc","nocol")

dat3[,vars_rev] <- dapply(dat3[,vars_rev], MARGIN = 2,
                    function(x) carrec(x,"0 = 'No'; 1 = 'Yes' "))

head(dat3)

# Inspect frequencies for categorical indicators.
freq(dat3[,1:4])

# Describe continuous achievement indicators.
summarytools::descr(dat3[,5:8])

# Fit mixed LCA models with continuous and categorical indicators.
fit3 <- lca(data = dat3, nclasses = 1:7,
           gaussian = c("ach9","ach10","ach11","ach12"),
           multinomial = c("hm","he","voc","nocol"))

# Extract and display fit indices for the mixed LCA models.
fit_ind3 <- getfit(fit3, digits = 2)
fit_ind3

# Plot model-fit indices for mixed LCA models.
plot(fit_ind3)

# Select the 2-class mixed model and inspect profile output.
fit3_2 <- fit3[[2]]

pr2 <- latInspect(fit3_2, "profile")
pr2$class
pr2$item[1:4]

# Inspect remaining profile output for the mixed model.
pr2$item[5:8]

# Plot mixed-model profiles.
plots_f3 <- plot(fit3_2)

plots_f3[[1]]/plots_f3[[2]]


# ---- Mixed LCA classification diagnostics and residuals ----
# Compute classification diagnostics for the mixed model.
cl_diag_mix <- lclass_diag(fit3_2, type = "all", digits = 3)

# Show posterior class totals.
cl_diag_mix$Sum.Posterior

# Show most-likely class totals.
cl_diag_mix$Sum.Mostlikely

# Show average posterior probabilities.
cl_diag_mix$AvePP

# Show entropy-based classification quality.
cl_diag_mix$Entropy

# Display classification diagnostics object.
cl_diag_cont

# Compute bivariate residual diagnostics for the mixed model.
res_diag3 <- lbvr(fit3_2, digits = 3)
res_diag3


# ---- Adding predictors to an LCA model ----
# Refit the 3-class categorical LCA model before adding predictors.
fit4 <- lca(data = dat, nclasses = 3,
           multinomial = c("contplt", "badge", "sgnptit", "bctprd",
                           "donprty", "pbldmna", "volunfp","pstplonl") )

# Add gender and age predictors to the 3-class model.
fit4_step2 <- lca(data = dat,
                  nclasses = 3,
                  multinomial = c("contplt", "badge",
                                  "sgnptit", "bctprd",
                                  "donprty", "pbldmna",
                                  "volunfp","pstplonl"),
                 X = c("gndr", "agea"),
                 model = fit4)

# Inspect predictor coefficients.
latInspect(fit4_step2, "coefs")

# Compute standard errors for predictor coefficients.
se_step2 <- se(fit4_step2)
se_step2$table$beta

# Compute 95% confidence intervals for predictor coefficients.
ci_step2 <- ci(fit4_step2, confidence = 0.95)
ci_step2$table$beta

# Plot the Odds Ratios + confidence intervals in forestplot style:
x <- plot_coeffs(fit4_step2,
                 what = "OR",        # Odds ratios
                 effects = "coding", # Effects-coding parametrization
                 confidence = 0.95,  # Confidence level
                 xlim = c(0, 3),     # Plot limits for x-axis
                 mfrow = c(2, 2))    # Grid of viewer panel

# Preview predicted class probabilities.
head(predict(fit4_step2))

# Predict class probabilities for gender profiles at age 50.
predict(fit4_step2,
        new = rbind(c(0,50),
                    c(1,50)))

# Predict class probabilities for selected gender and age profiles.
predict(fit4_step2,
        new = rbind(c(0,45),
                    c(0,50),
                    c(1,45),
                    c(1,50)))

control = list(opt     = "lbfgs", # newton
               maxit   = 1000L,   # Maximum number of iterations
               rstarts = 16L,     # Number of random starts
               cores   = 1L,      # Number or parallel runs
               eps     = 1e-05)   # Convergence criteria (gradient norm)

latInspect(fit4_step2, what = "convergence")

fit4_step2 <- lca(data = dat,
                  nclasses = 3,
                  multinomial = c("contplt", "badge",
                                  "sgnptit", "bctprd",
                                  "donprty", "pbldmna",
                                  "volunfp","pstplonl"),
                  X = c("gndr", "agea"),
                  model = fit4,
                  control = control)

latInspect(fit4_step2, what = "convergence")

penalty_defaults <- list(
  class = list(alpha = 1), # Penalty for class probabilities
  prob  = list(alpha = 1), # Penalty for item probabilities (multinomial items)
  var   = list(alpha = 1), # Penalty for item variances (gaussian items)
  Sigma = list(alpha = 1)  # Penalty for covariance matrices (gaussian items)
)

latInspect(fit4, what = "profile")
set.seed(15)
fit4.2 <- lca(data = dat, nclasses = 3,
              multinomial = c("contplt", "badge", "sgnptit", "bctprd",
                              "donprty", "pbldmna", "volunfp","pstplonl"),
              penalties = list(prob = list(alpha = 10)))
latInspect(fit4.2, what = "profile")


# Default arguments for the optimizer.
control = list(opt     = "lbfgs", # newton
               maxit   = 1000L,   # Maximum number of iterations
               rstarts = 16L,     # Number of random starts
               cores   = 1L,      # Number or parallel runs
               eps     = 1e-05)   # Convergence criteria (gradient norm)



