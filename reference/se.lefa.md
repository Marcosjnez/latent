# Standard Errors for Exploratory Factor Analysis

Standard Errors for Exploratory Factor Analysis

## Usage

``` r
# S3 method for class 'lefa'
se(x, type = NULL, parameters = NULL, digits = 4L, ..., fit = x)
```

## Arguments

- x:

  Fitted model or collection of fitted models passed to the method.

- type:

  Optional compatibility argument passed to
  [`se.multistep()`](https://marcosjnez.github.io/latent/reference/se.multistep.md).

- parameters:

  Optional parameter specification. By default, standard errors are
  returned for the rotated loadings, factor covariance matrices, and
  factor means.

- digits:

  Non-negative integer used to format parameter tables. If `NULL`,
  values are not rounded.

- ...:

  Additional arguments passed to
  [`se.multistep()`](https://marcosjnez.github.io/latent/reference/se.multistep.md).

- fit:

  A fitted object inheriting from class `"lefa"`.

## Value

The result of
[`se.multistep()`](https://marcosjnez.github.io/latent/reference/se.multistep.md).
