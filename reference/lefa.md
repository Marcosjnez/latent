# Exploratory Factor Analysis

Fit an exploratory factor analysis model by first estimating an
orthogonal factor model with
[`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md) and
subsequently rotating the fitted model with
[`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).

## Usage

``` r
lefa(data = NULL, nfactors = 1L, estimator = "ml",
     projection = "oblq", rotation = "oblimin",
     model = NULL, ordered = FALSE, group = NULL,
     sample.cov = NULL, sample.mean = NULL, sample.nobs = NULL,
     positive = FALSE, penalties = TRUE,
     missing = "pairwise.complete.obs",
     std.lv = TRUE, std.ov = TRUE,
     meanstructure = FALSE,
     parameterization = NULL,
     likelihood = NULL, se = TRUE,
     message = FALSE, do.fit = TRUE,
     mimic = "latent", control = NULL,
     rotation.control = NULL, ...)
```

## Arguments

- data:

  Optional data frame or matrix containing the observed variables.
  Alternatively, sample.cov can be supplied.

- nfactors:

  Integer. Number of factors used when `model = NULL`.

- estimator:

  Estimation method passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- projection:

  Rotation projection passed to
  [`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).
  Available options are `"orth"`, `"oblq"`, and `"poblq"`.

- rotation:

  Rotation criterion passed to
  [`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).

- model:

  Optional lavaan model syntax. If `NULL`, an exploratory lower-diagonal
  loading model is generated automatically.

- ordered:

  Logical value indicating whether indicators are ordinal. The character
  value `"yule"` requests Yule correlations.

- group:

  Optional character string identifying the grouping variable.

- sample.cov:

  Optional sample covariance matrix or list of covariance matrices
  passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- sample.mean:

  Optional sample mean vector or list of vectors passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- sample.nobs:

  Optional number of observations passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- positive:

  Logical. Request the positive-definite parameterization used by
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- penalties:

  Logical value or list controlling regularization in
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- missing:

  Missing-data method passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- std.lv:

  Logical. Standardize latent variables in the unrotated model.

- std.ov:

  Logical. Standardize observed variables in the unrotated model.

- meanstructure:

  Logical. Estimate the observed-variable mean structure.

- parameterization:

  Optional parameterization passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- likelihood:

  Character string controlling the likelihood convention in
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- se:

  Logical or character controlling standard errors in
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- message:

  Logical. Print progress messages during CFA estimation.

- do.fit:

  Logical. If `FALSE`, return the unrotated `lcfa` model specification
  without fitting or rotation.

- mimic:

  Retained for backward compatibility. Only `"latent"` is currently
  supported.

- control:

  Optional list of optimization controls passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- rotation.control:

  Optional list of optimization controls passed to
  [`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).

- ...:

  Additional arguments. CFA/lavaan arguments are passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md);
  arguments required by the selected rotation criterion or projection
  are passed only to
  [`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).

## Value

A fitted object of class `"lefa"`. The unrotated `lcfa` object is stored
in its `extra` slot. If `do.fit = FALSE`, the unfitted `lcfa`
specification is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lefa(data = HolzingerSwineford1939,
            nfactors = 3L,
            rotation = "oblimin")
} # }
```
