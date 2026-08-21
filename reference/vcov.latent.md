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

A list containing the selected variance-covariance matrix, standard
errors, and cumulative transformation Jacobian. Jacobian rows correspond
to the selected transformed parameters and columns to the freely
estimated parameter coordinates.
