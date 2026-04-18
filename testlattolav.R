

#### CFA ####

library(latent)
library(lavaan)

model <- 'visual  =~ x1 + x2 + x3
          textual =~ x4 + x5 + x6
          speed   =~ x7 + x8 + x9'

set.seed(2026)

#######
#######
####### continuous indicators
#######
#######

#####
## single group comparison
#####

fit <- lcfa(data=HolzingerSwineford1939, model = model,
            estimator = "ml", positive = FALSE, se = T,
            ordered = FALSE, acov = "standard",
            std.lv = TRUE, std.ov = FALSE,
            mimic = "latent", do.fit = TRUE,
            meanstructure = T, likelihood = "normal",
            control = NULL)

# With lavaan:
fit2 <- lavaan::cfa(model, data = HolzingerSwineford1939,
                    estimator = "ml",
                    likelihood = "normal",
                    std.lv = TRUE, std.ov = FALSE,
                    meanstructure = T)

fit2@loglik$loglik
fit@loglik

fit@loss
fit2@Fit@fx*2

latInspect(fit, what = "fit")


# transform lcfa object
fit_lav <- lcfa_to_lavaan(fit)

## comapre fit indices
round(cbind(fitMeasures(fit_lav),
            fitMeasures(fit2)),3)

## compare lavresiduals
lavResiduals(fit_lav)
lavResiduals(fit2)

## compare modindices
modindices(fit_lav, sort = T)[1:10,]
modindices(fit2, sort = T)[1:10,]

summary(fit_lav,
        standardized = T,
        rsquare = T,
        fit.measures= T)

summary(fit2,
        standardized = T,
        rsquare = T,
        fit.measures= T)



#####
## multiple group comparison
#####

fitmg <- lcfa(data=HolzingerSwineford1939, model = model, group = "school",
            estimator = "ml", positive = FALSE, se = T,
            ordered = FALSE, acov = "standard",
            std.lv = TRUE, std.ov = FALSE,
            mimic = "latent", do.fit = TRUE,
            meanstructure = T, likelihood = "normal",
            control = NULL)

# With lavaan:
fit2mg <- lavaan::cfa(model, data = HolzingerSwineford1939, group = "school",
                    estimator = "ml",
                    likelihood = "normal",
                    std.lv = TRUE, std.ov = FALSE,
                    meanstructure = T)

fit2mg@loglik$loglik
fitmg@loglik

fitmg@loss
fit2mg@Fit@fx*2

fitmg@Optim$f
latInspect(fitmg, what = "fit")


# transform lcfa object
fit_lavmg <- lcfa_to_lavaan(fitmg)

## comapre fit indices
round(cbind(fitMeasures(fit_lavmg),
            fitMeasures(fit2mg)),3)

## compare lavresiduals
lavResiduals(fit_lavmg)
lavResiduals(fit2mg)

## compare modindices
modindices(fit_lavmg, sort = T)[1:10,]
modindices(fit2mg, sort = T)[1:10,]

summary(fit_lavmg,
        standardized = T,
        rsquare = T,
        fit.measures= T)

summary(fit2mg,
        standardized = T,
        rsquare = T,
        fit.measures= T)



#####
## single group comparison missing data
#####

data_missing <- HolzingerSwineford1939
data_missing[1:19, "x1"] <- NA

fit_miss <- lcfa(data=data_missing, model = model,
                 estimator = "ml",
                 positive = FALSE,
                 ordered = FALSE,
                 std.lv = FALSE,
                 # missing = "pairwise.complete.obs",
                 missing = "fiml",
                 do.fit = TRUE,
                 std.ov = FALSE,
                 # meanstructure = TRUE,
                 control = NULL)

# With lavaan:
fit2_miss <- lavaan::cfa(model, data = data_missing,
                    estimator = "ml",
                    likelihood = "normal",
                    std.lv = TRUE, std.ov = FALSE,
                    meanstructure = T, missing = "fiml")

fit2_miss@loglik$loglik
fit_miss@loglik

fit_miss@loss
fit2_miss@Fit@fx

latInspect(fit_miss, what = "fit")


# transform lcfa object
fit_lav_miss <- lcfa_to_lavaan(fit_miss)

## comapre fit indices
round(cbind(fitMeasures(fit_lav_miss),
            fitMeasures(fit2_miss)),3)

## compare lavresiduals
lavResiduals(fit_lav_miss)
lavResiduals(fit2_miss)

## compare modindices
modindices(fit_lav_miss, sort = T)[1:10,]
modindices(fit2_miss, sort = T)[1:10,]

summary(fit_lav_miss,
        standardized = T,
        rsquare = T,
        fit.measures= T)

summary(fit2_miss,
        standardized = T,
        rsquare = T,
        fit.measures= T)




#######
#######
####### categorical indicators
#######
#######



samples <- unique(hexaco$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(hexaco$sample == x))
names(Ns) <- samples

# Subset the items pertaining to the HEXACO-100
selection <- 5:104
full <- hexaco[, selection]

mooc <- full[hexaco$sample == samples[2], ]
dim(mooc)
mooc$school <- rep(c("s1", "s2"), times = c(2000, 2286))

model.EM <- "FEA =~ hexemfea146 + hexemfea170 + hexemfea74 + hexemfea2
             ANX =~ hexemanx128 + hexemanx8 + hexemanx80 + hexemanx176
             DEP =~ hexemdep62 + hexemdep182 + hexemdep134 + hexemdep158
             SEN =~ hexemsen44 + hexemsen164 + hexemsen20 + hexemsen68"

fitcat <- lcfa(model = model.EM, data = mooc,
            ordered = TRUE, estimator = "dwls",
            #group = "school",
            meanstructure = T,std.lv = TRUE,
            do.fit = TRUE, control = NULL)


# With lavaan:
fit2cat <- lavaan::cfa(model = model.EM, data = mooc,
                    ordered = TRUE,
                    estimator = "dwls",
                    #group = "school",
                    # likelihood = "wishart",
                    std.lv = TRUE, std.ov = F,
                    meanstructure = T,
                    parameterization = "delta")


fitcat@loss
fit2cat@Fit@fx*2

# transform lcfa object
fit_lavcat <- lcfa_to_lavaan(fitcat)

## comapre fit indices
round(cbind(fitMeasures(fit_lavcat),
            fitMeasures(fit2cat)),3)

## compare lavresiduals
lavResiduals(fit_lavcat)
lavResiduals(fit2cat)

## compare modindices
modindices(fit_lavcat, sort = T)[1:10,]
modindices(fit2cat, sort = T)[1:10,]

