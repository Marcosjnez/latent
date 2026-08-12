# Standard Errors for Latent Class Models

`se.llca()` computes standard or robust standard errors for selected
parameters of a fitted latent class model.

Covariance matrices are obtained through
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md),
which also propagates uncertainty from earlier estimation stages when
the fitted model contains a nested `"latent"` object.

## Usage

``` r
# S3 method for class 'llca'
se(fit, type = "standard", parameters = NULL, digits = 4L, ...)
```

## Arguments

- fit:

  A fitted object of class `"llca"`.

- type:

  Character string specifying the covariance estimator. Available
  options are `"standard"` for Hessian-based standard errors and
  `"robust"` for a LatentGold-style sandwich estimator.

- parameters:

  Optional parameter specification identifying the parameters for which
  standard errors should be returned. Labels must occur in
  `fit@modelInfo$transparameters_labels`. If `NULL`, the parameter
  blocks corresponding to the fitted model parameters are used.

- digits:

  Non-negative integer specifying the number of decimal places used in
  the formatted parameter tables. If `NULL`, values are not rounded.

- ...:

  Additional arguments passed to other methods.

## Value

A list containing:

- `table`:

  Parameter estimates and standard errors arranged according to the
  model parameter blocks.

- `table_se`:

  Standard errors in the same parameter structure as the fitted model.

- `se`:

  Named numeric vector of standard errors.

- `vcov`:

  Variance-covariance matrix of the selected parameters.

- `B`:

  Additional uncertainty component for sequential models, or an empty
  matrix for ordinary one-step models.

- `H`:

  Hessian or information matrix used in the calculation.

- `newH`:

  Corrected joint precision matrix when applicable.

- `jacob`:

  Jacobian matrix used for covariance propagation to transformed
  parameters.

## Details

Compute standard errors, covariance matrices, Hessians, and
transformation Jacobians for fitted latent class models.

With `type = "standard"`, covariance matrices are based on the inverse
Hessian of the fitted objective function.

With `type = "robust"`, the empirical covariance of case or
response-pattern score contributions is combined with the Hessian to
produce a sandwich covariance estimator.

When a previous fitted `"latent"` object is stored in the model
specification of `fit`,
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md)
propagates uncertainty from that earlier estimation step. Nested
sequential models are processed recursively.

Standard errors for transformed parameters are obtained using the delta
method and the Jacobians of the corresponding parameter transformations.

## See also

[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md),
[`hessian.latent`](https://marcosjnez.github.io/latent/reference/hessian.latent.md),
[`jacobian.latent`](https://marcosjnez.github.io/latent/reference/jacobian.latent.md),
`ci`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lca(data = empathy, nclasses = 3L,
           gaussian = c("ec1", "ec2", "ec3"))

se(fit)
se(fit, type = "robust")

# Standard errors for selected transformed parameters
se(fit, parameters = fit@modelInfo$trans[c("class", "beta")])
} # }
```
