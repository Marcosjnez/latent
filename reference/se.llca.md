# Standard errors for latent class models

Computes standard errors and variance-covariance matrices for fitted
latent class models.

## Usage

``` r
# S3 method for class 'llca'
se(fit, type = "standard", digits = 3L, ...)
```

## Arguments

- fit:

  A fitted object of class `"llca"`.

- type:

  Character string indicating the standard-error estimator. Available
  options are `"standard"` for Hessian-based standard errors and
  `"robust"` for the LatentGold-style sandwich estimator.

- digits:

  Non-negative integer indicating the number of decimal places used in
  the formatted table. If `NULL`, the table is returned without
  rounding.

- ...:

  Additional arguments passed to other methods.

## Value

A list with the following components:

- `table`:

  A list of formatted parameter tables containing estimates and standard
  errors.

- `table_se`:

  A list containing the standard errors arranged in the same parameter
  structure as the fitted model.

- `se`:

  A named numeric vector of standard errors.

- `vcov`:

  The variance-covariance matrix of the model parameters.

- `B`:

  The empirical score covariance matrix. It is an empty matrix for
  ordinary Hessian-based standard errors and may be `NULL` when it is
  not applicable.

- `H`:

  The Hessian matrix, when available.

- `newH`:

  The adjusted Hessian used by the robust estimator, when available.

## Details

For a regular one-step model, `type = "standard"` obtains the covariance
matrix from the inverse Hessian. With `type = "robust"`, the covariance
matrix is computed from a sandwich estimator using the score
contribution of each observed response pattern.

When the fitted model contains a previous `"llca"` object in
`fit@modelInfo$control_optimizer$model`, the standard errors are
adjusted for two-step estimation through
[`se_twostep()`](https://marcosjnez.github.io/latent/reference/se_twostep.md).

The `digits` argument affects only the formatted `table`. The numeric
standard errors and covariance matrices are returned without rounding.

## See also

`ci()`, [`vcov()`](https://rdrr.io/r/stats/vcov.html)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lca(data = empathy, nclasses = 3L,
           gaussian = c("ec1", "ec2", "ec3"))

se(fit)
se(fit, type = "robust", digits = 4L)
} # }
```
