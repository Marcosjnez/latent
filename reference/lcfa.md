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
     do.fit = TRUE, control.moments = NULL, sort = TRUE, ...)
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

  Logical. Standardize observed variables. With direct or
  saturated-moment FIML, raw variables are standardized within
  substantive groups before missingness patterns are constructed;
  observed-variable means remain freely estimated.

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

  Optional list of CFA optimization controls. These are separate from
  `control.moments`.

- message:

  Logical. Print progress messages.

- do.fit:

  Logical. If `FALSE`, return the prepared but unfitted `"lcfa"` or
  `"multistep_lcfa"` object.

- control.moments:

  Optional named list passed as `control` to the sample-moment
  estimators
  ([`lpoly()`](https://marcosjnez.github.io/latent/reference/lpoly.md),
  [`lpearson()`](https://marcosjnez.github.io/latent/reference/lpearson.md),
  [`lmean()`](https://marcosjnez.github.io/latent/reference/lmean.md),
  [`lyule()`](https://marcosjnez.github.io/latent/reference/lyule.md),
  or
  [`lmvnorm()`](https://marcosjnez.github.io/latent/reference/lmvnorm.md)).
  For example, `list(cores = 4L)` requests four OpenMP threads for
  `polyfast()` polychoric estimation. The two-step ACOV does not use
  `cores` or explicit OpenMP. NULL uses the moment estimators' own
  defaults, independently of `control`. Group suffixes and uncertainty
  propagation are set internally. This argument has no effect when no
  separate moment estimator is fitted (for example, supplied sample
  moments or direct FIML).

- sort:

  Logical. Sort factors by decreasing variance-adjusted sums of squared
  loadings and make the largest absolute loading positive. The default
  is TRUE. If any modeled loading is fixed, including a loading fixed
  for factor-scale identification, sorting is automatically disabled so
  the fitted factor order and loading signs remain exactly as specified.
  Factor identities and the native estimation constraints are retained;
  see
  [`sort_factors()`](https://marcosjnez.github.io/latent/reference/sort_latent.md).
  FALSE leaves the output unchanged.

- ...:

  Additional arguments passed to lavaan and the sample-statistic
  estimators where applicable.

## Value

An S4 object of class `"lcfa"` for ordinary likelihood analyses, or
`"multistep_lcfa"` when at least one parent object is marked for
uncertainty propagation.

## Details

The model-implied observed means are computed as
\$\$\widehat{\mu}=\nu+\Lambda\alpha,\$\$ where \\\nu\\ contains
observed-variable intercepts and \\\alpha\\ contains latent-factor
means. For ordinal models, standardized model thresholds are computed
from the unstandardized thresholds, model-implied means, and
model-implied variances.

Direct FIML creates one likelihood contribution for every missingness
pattern and substantive group. Saturated-moment FIML instead stores one
`lmvnorm` source object in `extra`; its uncertainty is propagated
automatically.

## Examples

``` r
if (FALSE) { # \dontrun{
HS.model <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
'
fit <- lcfa(model = HS.model, data = HolzingerSwineford1939)
} # }
```
