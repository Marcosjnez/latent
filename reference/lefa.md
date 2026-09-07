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
     mimic = "latent", control.efa = NULL,
     control.rotation = NULL, control.moments = NULL, sort = TRUE, ...)
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

  A criterion name, character vector, or named list of component
  parameter lists passed to
  [`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).
  Multiple components are summed. Each list component can select `items`
  and `factors`; repeated criterion names are supported. See
  [`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md)
  for details.

- model:

  Optional lavaan model syntax. If `NULL`, an exploratory loading model
  is generated automatically. By default, the loading matrix is dense,
  its columns are mutually orthogonal, and their norms are unrestricted
  (`control.efa$orth.lambda = TRUE`). Set this control to FALSE to use a
  lower-triangular loading matrix. This option requires `std.lv = TRUE`
  and `positive = FALSE`.

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

  Logical. Standardize latent variables in the unrotated model. It must
  be `TRUE` when `control.efa$orth.lambda = TRUE`, so the factor
  covariance matrix is fixed to the identity and factor scale is
  represented by the unrestricted column norms of the loading matrix.

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

- control.efa:

  Optional list of controls passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md). The
  defaults are `rstarts = 3L`, `se_method = "KKT"`, and
  `orth.lambda = TRUE`.

- control.rotation:

  Optional list of controls passed to
  [`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).
  The defaults are `rstarts = 10L` and `se_method = "KKT"`.

- control.moments:

  Optional named list passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md)'s
  sample-moment estimators, independently of `control.efa` and
  `control.rotation`. For example, `list(cores = 4L)` controls
  `polyfast()` polychoric estimation. The two-step ACOV does not use
  `cores` or explicit OpenMP. NULL uses the moment estimators' defaults.
  See [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md)
  for details.

- sort:

  Logical. Sort the final rotated factors and choose their signs using
  [`sort_factors()`](https://marcosjnez.github.io/latent/reference/sort_latent.md).
  Defaults to TRUE. When any component is a target criterion or selects
  `factors`, factor order is retained and only signs are oriented. The
  intermediate unrotated model retains its fitting coordinates.

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

## Details

Two equivalent identification schemes are available for the unrotated
EFA model. With `control.efa$orth.lambda = FALSE`, the model uses a
lower-triangular loading matrix and an identity factor covariance
matrix. The default, `control.efa$orth.lambda = TRUE`, makes every
loading free. The factor covariance matrix remains the identity, and
only the off-diagonal elements of \\\Lambda^\top\Lambda\\ are
constrained to zero. The column norms of \\\Lambda\\ remain free. The
\\q(q-1)/2\\ orthogonality constraints replace the same number of
lower-triangular zeros, so both schemes have the same effective number
of parameters.

With `projection = "poblq"`, either `constraints` or `oblique` must be
supplied through `...`. The former uses arbitrary structural
constraints, whereas the latter gives the sizes of consecutive oblique
blocks; any remaining factors form one orthogonal block. They cannot be
used together.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lefa(data = HolzingerSwineford1939,
            nfactors = 3L,
            rotation = "oblimin")

fit_orth_lambda <- lefa(
  data = HolzingerSwineford1939,
  nfactors = 3L,
  rotation = "oblimin",
  control.efa = list(orth.lambda = TRUE)
)
} # }
```
