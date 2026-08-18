# Confirmatory Factor Analysis

Fit confirmatory factor analysis models using lavaan model syntax and
the optimization infrastructure of latent.

## Usage

``` r
lcfa(data = NULL, model = NULL, estimator = "ml",
     ordered = FALSE, group = NULL,
     sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
     positive = FALSE, penalties = FALSE,
     missing = "pairwise.complete.obs",
     std.lv = FALSE, std.ov = FALSE,
     acov = "standard", meanstructure = TRUE,
     parameterization = NULL,
     likelihood = NULL, se = TRUE,
     control = NULL, message = FALSE,
     do.fit = TRUE, ...)
```

## Arguments

- data:

  Optional data frame or matrix containing the observed variables. If
  NULL, sample.cov and sample.nobs must be supplied.

- model:

  Confirmatory factor model specified using lavaan syntax.

- estimator:

  Estimation method. Available options include `"ml"`, `"uls"`, and
  `"dwls"`.

- ordered:

  Logical value indicating whether indicators are ordinal. The character
  value `"yule"` requests Yule correlations.

- group:

  Optional character string identifying the grouping variable.

- sample.cov:

  Optional sample covariance matrix or list of covariance matrices. Used
  when data is NULL.

- sample.mean:

  Optional sample mean vector or list of vectors. Required when data is
  NULL and meanstructure = TRUE.

- sample.nobs:

  Optional number of observations, or one value per group, used when
  data is NULL.

- positive:

  Logical. If `TRUE`, positive-definite covariance structures are
  imposed through the corresponding manifold parameterization.

- penalties:

  Logical value or list controlling regularization.

- missing:

  Missing-data method.

- std.lv:

  Logical. Standardize latent variables.

- std.ov:

  Logical. Standardize observed variables.

- acov:

  Character string selecting the variance-covariance estimator for the
  sample statistics. Available options are `"standard"` and `"robust"`.

- meanstructure:

  Logical. Estimate the observed-variable mean structure.

- parameterization:

  Optional parameterization specification.

- likelihood:

  Character string controlling the normal/Wishart likelihood convention.

- se:

  Logical. Compute standard errors before returning the fitted object.

- control:

  Optional list of optimization controls.

- message:

  Logical. Print progress messages.

- do.fit:

  Logical. If `FALSE`, return the prepared but unfitted `"lcfa"` object.

- ...:

  Additional arguments passed to lavaan and the sample-statistic
  estimators where applicable.

## Value

An S4 object of class `"lcfa"`, which inherits from `"latent"`.

## Examples

``` r
if (FALSE) { # \dontrun{
HS.model <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
'

fit <- lcfa(model = HS.model, data = HolzingerSwineford1939)
summary(fit, digits = 3L)

S <- cov(HolzingerSwineford1939[, paste0("x", 1:9)])
M <- colMeans(HolzingerSwineford1939[, paste0("x", 1:9)])
fit_cov <- lcfa(model = HS.model, sample.cov = S, sample.mean = M,
                sample.nobs = nrow(HolzingerSwineford1939))
} # }
```
