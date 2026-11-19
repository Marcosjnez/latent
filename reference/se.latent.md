# Standard Errors for Latent Objects

Standard Errors for Latent Objects

## Usage

``` r
# S3 method for class 'latent'
se(x, type = "information", parameters = NULL, digits = 4L, ..., fit = x)
```

## Arguments

- x:

  Fitted model or collection of fitted models passed to the method.

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

- fit:

  A fitted object inheriting from class `"latent"`.

## Value

A list containing parameter tables, standard errors, covariance
matrices, and derivative information.
