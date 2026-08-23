# Standard Errors for Exploratory Factor Analysis

Standard Errors for Exploratory Factor Analysis

## Usage

``` r
# S3 method for class 'lefa'
se(fit, type = NULL, parameters = NULL, digits = 4L, ...)
```

## Arguments

- fit:

  A fitted object inheriting from class `"lefa"`.

- type:

  Optional compatibility argument.

- parameters:

  Optional parameter specification. By default, standard errors are
  returned for the rotated loadings, factor covariance matrices, and
  factor means.

- digits:

  Non-negative integer used to format parameter tables. If `NULL`,
  values are not rounded.

- ...:

  Additional arguments reserved for future methods.

## Value

Standard errors propagated from the fitted unrotated `lcfa` model to the
rotated EFA parameters.
