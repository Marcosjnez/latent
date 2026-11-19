# Fit Indices for Exploratory Factor Analysis

Fit Indices for Exploratory Factor Analysis

## Usage

``` r
# S3 method for class 'lefa'
getfit(x, digits = 3L, fit_matrix = NULL, ..., model = x)
```

## Arguments

- x:

  Fitted model or collection of fitted models passed to the method.

- digits:

  Number of decimal places used in the returned values.

- fit_matrix:

  Optional precomputed internal fit matrix passed to
  [`getfit.lcfa()`](https://marcosjnez.github.io/latent/reference/getfit.lcfa.md).

- ...:

  Additional arguments to the generic.

- model:

  A fitted object inheriting from class `"lefa"`.

## Value

A named numeric vector containing model dimensions and fit statistics
from the unrotated model. Rotation does not change model fit.
