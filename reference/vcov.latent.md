# Variance-Covariance Matrix for Latent Objects

Variance-Covariance Matrix for Latent Objects

## Usage

``` r
# S3 method for class 'latent'
vcov(object, v, parameters = NULL, ..., fit = object)
```

## Arguments

- object:

  Fitted model or collection of fitted models passed to the method.

- v:

  Variance-covariance matrix of the freely estimated, untransformed
  parameters, in the original fitting coordinates.

- parameters:

  Optional parameter specification identifying the parameters or
  transformed parameters to return. The default follows the original
  order of the freely estimated parameters.

- ...:

  Additional arguments to the generic.

- fit:

  A fitted object inheriting from class `"latent"`.

## Value

A list containing the selected variance-covariance matrix and standard
errors.

## Details

Transformation covariances are propagated incrementally in C++. The
internal Jacobian workspace used by that computation is not returned
because it is not the cumulative Jacobian exposed by
[`jacobian.latent()`](https://marcosjnez.github.io/latent/reference/jacobian.latent.md).
