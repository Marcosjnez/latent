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
     control.rotation = NULL, control.moments = NULL, ...)
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
  Available options are `"orth"`, `"oblq"`, and `"poblq"`. With `"oblq"`
  (the default), a non-NULL `oblique` argument supplied through `...`
  selects `"poblq"` internally. If `oblique` is omitted or NULL, all
  factors are allowed to be oblique.

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

- ...:

  Additional arguments. CFA/lavaan arguments are passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md);
  arguments required by the selected rotation criterion or projection
  are passed only to
  [`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md).
  This includes the optional `oblique` block specification described in
  Details, which is never passed to
  [`lcfa()`](https://marcosjnez.github.io/latent/reference/lcfa.md).

## Value

A fitted object of class `"lefa"`. The unrotated `lcfa` object is stored
in its `extra` slot. If `do.fit = FALSE`, the unfitted `lcfa`
specification is returned.

## Details

Fitted objects retain the estimated factor order and signs. Sorting is
available only through `latInspect(fit, sort = TRUE)` and does not
modify the fitted object or its standard-error calculations.

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

**The oblique argument.** Supply `oblique` through `...`. With
`projection = "oblq"` (the default), an omitted or NULL `oblique`
retains the fully oblique projection. A non-NULL block specification
selects `"poblq"` internally before CFA and rotation arguments are
separated. The resulting rotation is the same as explicitly requesting
`projection = "poblq"` with that block specification.

A numeric vector specifies the sizes of consecutive oblique blocks; each
size must be a positive integer and their sum cannot exceed the number
of factors. For example, `oblique = c(2, 2)` defines blocks 1–2 and 3–4.
A list instead gives explicit one-based factor positions, such as
`oblique = list(c(1, 3), c(2, 4))`. Blocks must be non-empty and
disjoint, and positions must lie within the factor range; singleton
blocks are allowed. Any unlisted factors are mutually orthogonal and
orthogonal to every listed block. Factor order and target/weight
indexing are unchanged. Invalid block specifications are rejected when
the rotation model is constructed, rather than reverting to a fully
oblique rotation.

Numeric block-size vectors and flat lists of factor-position vectors are
shared by all groups. For group-specific blocks, use an outer list with
one block list per group, in input group order; for example,
`list(list(c(1, 3), c(2, 4)), list(1:2, 3:4))` for two groups. See
[`lrotate()`](https://marcosjnez.github.io/latent/reference/lrotate.md)
for further details about group-specific arguments.

Explicit `projection = "poblq"` still requires either `oblique` or a
`constraints` matrix supplied through `...`. To use a `constraints`
matrix for rotation, select `"poblq"` explicitly. These two arguments
cannot be used together. A non-NULL `oblique` is not allowed with
`projection = "orth"`. As before, `do.fit = FALSE` returns only the
unrotated CFA specification; no rotation model is constructed in that
case.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lefa(data = HolzingerSwineford1939,
            nfactors = 3L,
            rotation = "oblimin")

# Factors 1 and 2 are oblique; factor 3 is orthogonal to both.
# projection = "poblq" is selected internally.
block_fit <- lefa(data = HolzingerSwineford1939,
                  nfactors = 3L, oblique = 2)

# Nonconsecutive positions: factors 1 and 3 may correlate.
position_fit <- lefa(data = HolzingerSwineford1939,
                     nfactors = 3L, projection = "oblq",
                     oblique = list(c(1, 3)))

fit_orth_lambda <- lefa(
  data = HolzingerSwineford1939,
  nfactors = 3L,
  rotation = "oblimin",
  control.efa = list(orth.lambda = TRUE)
)
} # }
```
