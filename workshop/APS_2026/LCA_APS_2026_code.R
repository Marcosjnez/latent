# Latent Class Analysis: latent R package
# Workshop code

# Please, reinstall the last version of latent before starting:
# devtools::install_github("marcosjnez/latent", force = TRUE)

# Visit https://github.com/Marcosjnez/latent for more installation tips

# ================================================================
# Measurement models
# ================================================================
# Load packages, create scale scores, fit the CFA model, and plot the path diagram.

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

# Compare correlations among average item scores.

## average items
round(cor(dat[,c("vis","tex","spd")]),3)

# Inspect correlations among latent factors from the fitted CFA model.

## measurement model
lavInspect(fit, "cor.lv")

# ================================================================
# The latent package
# ================================================================
# Load the package used for the latent class analysis examples.

library(latent)

# ================================================================
# LCA with categorical indicators
# ================================================================
# Load helper packages, select the French ESS sample, recode indicators, and inspect frequencies.

library(collapse)
library(summarytools)
library(stevemisc)
library(ggplot2)
library(patchwork)

dat <- sbt(ESS_datLCA, cntry == "FR")
dat <- tidyr::drop_na(dat)

dat[,1:8] <- dapply(dat[,1:8], MARGIN = 2,
                    function(x) carrec(x,"1='Yes';2='No'"))

dat$gndr <- carrec(dat$gndr, "1='Men'; 2 ='Women' ")

freq(dat[,1:2])

# Inspect frequencies for the next pair of categorical indicators.

freq(dat[,3:4])

# Inspect frequencies for the following pair of categorical indicators.

freq(dat[,5:6])

# Inspect frequencies for the final pair of categorical indicators.

freq(dat[,7:8])

# Estimate categorical LCA models with one to seven classes.

fit1 <- lca(data = dat, nclasses = 1:7,
            multinomial = c("contplt", "badge", "sgnptit", "bctprd",
                            "donprty", "pbldmna", "volunfp","pstplonl") )

# Extract and inspect model-fit indices.

fit_ind <- getfit(fit1, digits = 2)
fit_ind

# Plot fit indices across the class solutions.

plot(fit_ind)

# Inspect the class proportions and first item probabilities for the 2-class solution.

pr2 <- latInspect(fit1[[2]], "profile")
pr2$class
pr2$item[1:3]

# Inspect the remaining item probabilities for the 2-class solution.

pr2$item[4:8]

# Inspect the class proportions and first item probabilities for the 3-class solution.

pr3 <- latInspect(fit1[[3]], "profile")
pr3$class
pr3$item[1:3]

# Inspect the remaining item probabilities for the 3-class solution.

pr3$item[4:8]

# Select and plot the 3-class solution.

fit1_3 <- fit1[[3]]

plot(fit1_3)

# Inspect estimated class proportions for the selected solution.

latInspect(fit1_3, "classes")

# ================================================================
# Classification diagnostics for categorical LCA
# ================================================================
# Compute classification-quality diagnostics for the selected model.

cl_diag <- lclass_diag(fit1_3, digits = 3, type = "all")

# Sum posterior probabilities by class.

cl_diag$Sum.Posterior

# Count most-likely class assignments.

cl_diag$Sum.Mostlikely

# Average posterior probabilities for the assigned classes.

cl_diag$AvePP

# Entropy-based classification quality.

cl_diag$Entropy

# Print the full classification diagnostics object.

cl_diag

# ================================================================
# Predicted classes for categorical LCA
# ================================================================
# Inspect most-likely class assignments.

latInspect(fit1_3, "state")

# Add predicted classes to the data and cross-tabulate them by gender.

dat$State <- latInspect(fit1_3, "state")
ctable(dat$gndr, dat$State, prop="c")
ctable(dat$gndr, dat$State, prop="r")

# ================================================================
# Global model fit
# ================================================================
# Inspect likelihood-ratio fit statistics for the selected model.

getfit(fit1_3)[c("L2","dof","pvalue")]

# ================================================================
# Residual evaluation for categorical LCA
# ================================================================
# Compute bivariate residual diagnostics.

res_diag <- lbvr(fit1_3, digits = 3)
res_diag

# Refit the categorical LCA model with a residual dependency between two indicators.

fit1_3_resid <- lca(data = dat, nclasses = 3,
                    multinomial = c("contplt", "badge", "sgnptit", "bctprd",
                                    "donprty", "pbldmna", "volunfp","pstplonl"),
                    model = list("contplt ~~ volunfp"))

fit1_3_resid # Print the model

# ================================================================
# LCA with continuous indicators
# ================================================================
# Prepare Big Five data, reverse-code selected items, and compute scale scores.

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

# Estimate continuous-indicator LCA models with one to seven classes.

set.seed(1987)

fit2 <- lca(data = bfi, nclasses = 1:7,
           gaussian = c("Agreeableness", "Conscientiousness",
                        "Extraversion", "Neuroticism", "Openness" ) )

# Extract and inspect fit indices for the continuous LCA models.

fit_ind_cont <- getfit(fit2, digits = 2)
fit_ind_cont

# Plot fit indices across continuous class solutions.

plot(fit_ind_cont)

# Compare class proportions for candidate continuous models.

latInspect(fit2[[3]], "classes")
latInspect(fit2[[4]], "classes")
latInspect(fit2[[5]], "classes")

# Select the 4-class solution and inspect the first profile outputs.

fit2_4 <- fit2[[4]]

pr4 <- latInspect(fit2_4, "profile")
pr4$class
pr4$item[1:2]

# Inspect the remaining profile outputs for the 4-class solution.

pr4$item[3:5]

# Plot profiles for the selected continuous model.

plot(fit2_4)

# Inspect class proportions for the selected continuous model.

latInspect(fit2_4, "classes")

# ================================================================
# Classification diagnostics for continuous LCA
# ================================================================
# Compute classification-quality diagnostics for the continuous model.

cl_diag_cont <- lclass_diag(fit2_4, type = "all", digits = 3)

# Sum posterior probabilities by class.

cl_diag_cont$Sum.Posterior

# Count most-likely class assignments.

cl_diag_cont$Sum.Mostlikely

# Average posterior probabilities for the assigned classes.

cl_diag_cont$AvePP

# Entropy-based classification quality.

cl_diag_cont$Entropy

# Print the full classification diagnostics object.

cl_diag_cont

# ================================================================
# Predicted classes for continuous LCA
# ================================================================
# Inspect most-likely class assignments.

latInspect(fit2_4, "state")

# Add predicted classes to the data and cross-tabulate them by gender.

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

# ================================================================
# Residual evaluation for continuous LCA
# ================================================================
# Compute bivariate residual diagnostics.

res_diag2 <- lbvr(fit2_4, digits = 3)
res_diag2

# Refit the continuous model with a residual covariance between two indicators.

fit2_4_resid <- lca(data = bfi, nclasses = 4,
            gaussian = c("Agreeableness", "Conscientiousness",
                         "Extraversion", "Neuroticism", "Openness" ),
            model = list("Agreeableness ~~ Extraversion"))

getfit(fit2_4_resid)

# ================================================================
# LCA with continuous and categorical indicators
# ================================================================
# Prepare the mixed-indicator data and recode categorical variables.

dat3 <- hs_academic

vars_rev <- c("hm","he","voc","nocol")

dat3[,vars_rev] <- dapply(dat3[,vars_rev], MARGIN = 2,
                    function(x) carrec(x,"0 = 'No'; 1 = 'Yes' "))

head(dat3)

# Inspect frequencies for the categorical indicators.

freq(dat3[,1:4])

# Describe the continuous achievement indicators.

summarytools::descr(dat3[,5:8])

# Estimate mixed LCA models with one to seven classes.

fit3 <- lca(data = dat3, nclasses = 1:7,
           gaussian = c("ach9","ach10","ach11","ach12"),
           multinomial = c("hm","he","voc","nocol"))

# Extract and inspect model-fit indices for the mixed LCA models.

fit_ind3 <- getfit(fit3, digits = 2)
fit_ind3

# Plot fit indices across mixed class solutions.

plot(fit_ind3)

# Select the 2-class mixed model and inspect the first profile outputs.

fit3_2 <- fit3[[2]]

pr2 <- latInspect(fit3_2, "profile")
pr2$class
pr2$item[1:4]

# Inspect the remaining profile outputs for the mixed model.

pr2$item[5:8]

# Plot the categorical and continuous parts of the mixed profile.

plots_f3 <- plot(fit3_2)

plots_f3[[1]]/plots_f3[[2]]

# ================================================================
# Classification diagnostics for mixed LCA
# ================================================================
# Compute classification-quality diagnostics for the mixed model.

cl_diag_mix <- lclass_diag(fit3_2, type = "all", digits = 3)

# Sum posterior probabilities by class.

cl_diag_mix$Sum.Posterior

# Count most-likely class assignments.

cl_diag_mix$Sum.Mostlikely

# Average posterior probabilities for the assigned classes.

cl_diag_mix$AvePP

# Entropy-based classification quality.

cl_diag_mix$Entropy

# Print the full mixed-model classification diagnostics object.

cl_diag_mix

# ================================================================
# Residual evaluation for mixed LCA
# ================================================================
# Compute bivariate residual diagnostics for the mixed model.

res_diag3 <- lbvr(fit3_2, digits = 3)
res_diag3

# ================================================================
# Adding predictors to an LCA model
# ================================================================
# Refit the three-class categorical model before adding covariates.

fit4 <- lca(data = dat, nclasses = 3,
           multinomial = c("contplt", "badge", "sgnptit", "bctprd",
                           "donprty", "pbldmna", "volunfp","pstplonl") )

# Add gender and age as predictors of class membership.

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

# Compute standard errors and confidence intervals for predictor coefficients.

se_step2 <- se(fit4_step2)
se_step2$table$beta

# Compute 95% confidence intervals for predictor coefficients.
ci_step2 <- ci(fit4_step2, confidence = 0.95)
ci_step2$table$beta

# Plot odds ratios and confidence intervals for the predictor effects.

x <- plot_coeffs(fit4_step2,
                 what = "OR",        # Odds ratios
                 effects = "coding", # Effects-coding parametrization
                 confidence = 0.95,  # Confidence level
                 xlim = c(0, 3),     # Plot limits for x-axis
                 mfrow = c(2, 2))    # Grid of viewer panel

# Preview predicted class probabilities for the observed data.

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

# ================================================================
# Troubleshooting convergence
# ================================================================
# Fit a deliberately limited model and inspect log-likelihood, convergence, and gradients.

set.seed(2026)
fit <- lca(data = empathy,
           nclasses = 2L,
           gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
           control = list(maxit = 20))

latInspect(fit, what = "loglik")
latInspect(fit, what = "convergence")
latInspect(fit, what = "gradient")

# Fit a mixed-indicator model with one random start and inspect convergence diagnostics.

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

# Fit the same mixed-indicator model using a more complete optimizer-control specification.

fit <- lca(data = cancer,
           nclasses = 3L,
           gaussian = c("Age", "WeightIndex", "SystolicBloodPressure",
                        "DiastolicBloodPressure"),
           multinomial = c("PerformanceRating", "CardiovascularDiseaseHistory"),
           control = list(opt     = "lbfgs", # newton
                          maxit   = 1000L,   # Maximum number of iterations
                          rstarts = 16L,     # Number of random starts
                          cores   = 1L,      # Number or parallel runs
                          eps     = 1e-05)   # Convergence criteria (gradient norm))
)

# ================================================================
# Penalized estimation
# ================================================================
# Fit the empathy model without penalties to illustrate a problematic solution.

set.seed(59) # Pick a seed that creates the problem
fit_empathy <- lca(data = empathy,
           nclasses = 4L,
           gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
           penalties = FALSE) # Deactivate penalized-ML, just ML

latInspect(fit_empathy, what = "loglik")
latInspect(fit_empathy, what = "convergence")
latInspect(fit_empathy, what = "profile") # What's wrong in this profile?

# Refit the empathy model with the default penalized-ML settings.

set.seed(59) # Pick a seed that creates the problem
fit_empathy <- lca(data = empathy,
           nclasses = 4L,
           gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
           penalties = TRUE)

latInspect(fit_empathy, what = "loglik")
latInspect(fit_empathy, what = "convergence")
latInspect(fit_empathy, what = "profile")

# Store the default penalty settings used by the examples.

penalty_defaults <- list(
  class = list(alpha = 1), # Penalty for class probabilities
  prob  = list(alpha = 1), # Penalty for item probabilities (multinomial items)
  var   = list(alpha = 1), # Penalty for item variances (gaussian items)
  Sigma = list(alpha = 1)  # Penalty for covariance matrices (gaussian items)
)

# Apply a stronger variance penalty and compare the resulting profiles.

latInspect(fit_empathy, what = "profile") # Before
fit_empathy <- lca(data = empathy,
                   nclasses = 4L,
                   gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
                   penalties = list(var = list(alpha = 2))) # Extreme penalty
latInspect(fit_empathy, what = "profile") # After
