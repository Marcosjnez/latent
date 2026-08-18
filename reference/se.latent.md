# Standard Errors for Latent Objects

Standard Errors for Latent Objects

## Usage

``` r
# S3 method for class 'latent'
se(fit, type = "information", parameters = NULL, digits = 4L, ...)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

- type:

  Character string selecting `"information"` or `"robust"` covariance
  estimation.

- parameters:

  Optional parameter specification.

- digits:

  Non-negative integer used to format parameter tables, or `NULL` to
  avoid rounding.

- ...:

  Additional arguments.

## Value

A list containing parameter tables, standard errors, covariance
matrices, and derivative information.
