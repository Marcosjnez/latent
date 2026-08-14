# Exploratory Factor Analysis

Fit an exploratory factor analysis model by first estimating an
orthogonal factor model with
[`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md) and
subsequently rotating the estimated loading matrices with
[`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).

## Usage

``` r
lefa(data, nfactors = 1L, estimator = "ml",
     projection = "oblq", rotation = "oblimin",
     model = NULL, ordered = FALSE, group = NULL,
     sample.cov = NULL, nobs = NULL,
     positive = FALSE, penalties = TRUE,
     missing = "pairwise.complete.obs",
     std.lv = TRUE, do.fit = TRUE,
     mimic = "latent", control = NULL,
     ...)
```

## Arguments

- data:

  A data frame or matrix containing the observed variables.

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

- nobs:

  Optional number of observations passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- positive:

  Logical. Request the positive-definite parameterization used by
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- penalties:

  Logical value or list controlling regularization.

- missing:

  Missing-data method passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

- std.lv:

  Logical. Standardize latent variables in the unrotated model.

- do.fit:

  Logical. If `FALSE`, return the unrotated `lcfa` model specification
  without fitting or rotation.

- mimic:

  Retained for backward compatibility. Only `"latent"` is currently
  supported.

- control:

  Optional list of optimization controls passed to the unrotated
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md)
  model.

- ...:

  Additional arguments. CFA/lavaan arguments are passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md);
  arguments required by the selected rotation criterion or projection
  are passed only to
  [`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).

## Value

If the fitted model contains more than one factor, a list with
components `efa` (the unrotated `lcfa` fit) and `rotation` (the rotated
`latent` fit). For a one-factor model, the unrotated `lcfa` fit is
returned directly. If `do.fit = FALSE`, the unfitted `lcfa`
specification is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lefa(data = HolzingerSwineford1939,
            nfactors = 3L,
            rotation = "oblimin")
} # }
```
