# Variance-Covariance Matrix for Latent Objects

Variance-Covariance Matrix for Latent Objects

## Usage

``` r
# S3 method for class 'latent'
vcov(fit, v, parameters = NULL)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

- v:

  Variance-covariance matrix of the freely estimated, untransformed
  parameters.

- parameters:

  Optional parameter specification identifying the parameters or
  transformed parameters to return.

## Value

A list containing the selected variance-covariance matrix and standard
errors.

## Details

Transformation covariances are propagated incrementally in C++. The
internal Jacobian workspace used by that computation is not returned
because it is not the cumulative Jacobian exposed by
[`jacobian.latent()`](https://marcosjnez.github.io/latent/reference/jacobian.latent.md).
