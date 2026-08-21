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
     meanstructure = TRUE,
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
  `"dwls"`. The value `"fiml"` requests direct pattern-likelihood FIML
  unless `missing = "fiml"` is supplied explicitly.

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

  Missing-data method. `"ml"` uses direct pattern-likelihood FIML.
  `"fiml"` first estimates saturated incomplete-data moments with
  [`lmvnorm()`](https://marcosjnez.github.io/latent/reference/lmvnorm.md)
  and then fits CFA as a deterministic multistep estimator.

- std.lv:

  Logical. Standardize latent variables.

- std.ov:

  Logical. Standardize observed variables.

- meanstructure:

  Logical. Estimate the observed-variable mean structure.

- parameterization:

  Optional parameterization specification.

- likelihood:

  Character string controlling the normal/Wishart likelihood convention.

- se:

  Logical or character. `TRUE`, `"standard"`, and `"information"` use
  standard sampling covariance matrices for the sample statistics.
  `"robust"` requests robust sampling covariance matrices where
  implemented. `FALSE` skips computation of the final CFA standard
  errors. Ordinary ML and direct FIML use information from the
  likelihood; multistep analyses propagate the covariance of their
  parent statistics with
  [`se.multistep()`](https://marcosjnez.github.io/latent/reference/se.multistep.md).

- control:

  Optional list of optimization controls.

- message:

  Logical. Print progress messages.

- do.fit:

  Logical. If `FALSE`, return the prepared but unfitted
  `"multistep_lcfa"` object.

- ...:

  Additional arguments passed to lavaan and the sample-statistic
  estimators where applicable.

## Value

An S4 object of class `"lcfa"` for ordinary likelihood analyses, or
`"multistep_lcfa"` when at least one parent object is marked for
uncertainty propagation.

## Details

Direct FIML creates one likelihood contribution for every missingness
pattern and substantive group. Pattern means and covariance matrices are
fixed transformed parameters and are not stored as separate latent
objects. Saturated-moment FIML instead stores one `lmvnorm` source
object in `extra`; its uncertainty is propagated automatically.

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
