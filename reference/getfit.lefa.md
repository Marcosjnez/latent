# Fit Indices for Exploratory Factor Analysis

Fit Indices for Exploratory Factor Analysis

## Usage

``` r
# S3 method for class 'lefa'
getfit(model, digits = 3L, fit_matrix = NULL)
```

## Arguments

- model:

  A fitted object inheriting from class `"lefa"`.

- digits:

  Number of decimal places used in the returned values.

- fit_matrix:

  Optional precomputed internal fit matrix passed to
  [`getfit.lcfa()`](https://marcosjnez.github.io/latent/reference/getfit.lcfa.md).

## Value

A named numeric vector containing model dimensions and fit statistics
from the unrotated model. Rotation does not change model fit.
