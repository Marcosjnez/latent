# Jacobian Matrix for Latent Models

`jacobian.latent()` evaluates the derivatives of transformed parameters
at the fitted parameter estimates. The Jacobian is used by latent to
propagate covariance matrices through parameter transformations using
the delta method.

## Usage

``` r
# S3 method for class 'latent'
jacobian(fit, parameters = NULL)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

- parameters:

  Optional parameter specification identifying the transformed
  parameters for which derivatives should be returned. Parameter labels
  must occur in `fit@modelInfo$transparameters_labels`. If `NULL`, the
  parameter blocks corresponding to the fitted model parameters are
  used.

## Value

A numeric matrix containing the Jacobian for the selected parameters.
Row and column names correspond to transformed-parameter labels.

## Details

Compute the Jacobian matrix associated with parameter transformations in
a fitted latent variable model.

Only transformations required to obtain the selected parameters are
evaluated. Dependencies between transformations are identified
recursively, so parameters obtained through several successive
transformations are supported.

The same Jacobian machinery is used by
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md)
to apply the delta method to transformed parameters.

## See also

[`jacobian`](https://marcosjnez.github.io/latent/reference/jacobian.md),
[`hessian.latent`](https://marcosjnez.github.io/latent/reference/hessian.latent.md),
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md),
[`constraints_derivs.latent`](https://marcosjnez.github.io/latent/reference/constraints_derivs.latent.md)
