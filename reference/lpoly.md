# Polychoric Correlation Matrix

Estimate a polychoric correlation matrix from ordinal data using
one-step or two-step estimation.

## Usage

``` r
lpoly(data, method = "two-step", model = NULL,
      positive = FALSE, penalties = FALSE,
      start = NULL, do.fit = TRUE, message = FALSE,
      control = NULL, ...)
```

## Arguments

- data:

  A data frame or matrix containing ordinal variables coded numerically.

- method:

  Character string indicating the estimation method. Available options
  are `"one-step"` and `"two-step"`.

- model:

  Optional named list used with `method = "one-step"` to fix parameters
  or impose equality constraints. Names should correspond to parameter
  blocks such as `S` or the threshold blocks. Numeric values fix
  parameters, repeated character labels impose equality constraints, and
  `NA` leaves the corresponding parameter unchanged.

- positive:

  Logical. If `TRUE`, the correlation matrix is parameterized through an
  oblique manifold so that the resulting matrix is positive
  semidefinite. This option is available for one-step estimation.

- penalties:

  Logical value or named list controlling regularization. If `TRUE`, the
  default log-determinant penalty is used. Two-step estimation does not
  use penalties.

- start:

  Optional named list of starting values used with
  `method = "one-step"`. Names should correspond to parameter blocks in
  the model. Partial matrices/vectors and `NA` values are handled in the
  same way as in
  [`lca()`](https://marcosjnez.github.io/latent/reference/lca.md).

- do.fit:

  Logical. If `FALSE`, return the prepared but unfitted `"latent"`
  object.

- message:

  Logical. Print progress messages during estimation.

- control:

  Optional list of optimization controls.

- ...:

  Additional arguments reserved for future extensions.

## Value

An S4 object of class `"latent"`. The object contains the processed data
in `dataList`, the model and optimization structures in `modelInfo`,
optimizer output in `Optim`, and the estimated parameter structures in
`parameters` and `transformed_pars`.

## Details

With `method = "two-step"`, thresholds are obtained from the marginal
distributions and the polychoric correlations are then estimated
pairwise. With `method = "one-step"`, thresholds and correlations are
optimized jointly using the general optimization infrastructure of
latent.

The two-step method always uses the unconstrained pairwise estimator and
therefore does not support `model`, `start`, `positive`, or `penalties`.
Model constraints and custom starts are available for
`method = "one-step"`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lpoly(data = values)
fit_one_step <- lpoly(data = values, method = "one-step")
fit_positive <- lpoly(data = values, method = "one-step", positive = TRUE)
} # }
```
