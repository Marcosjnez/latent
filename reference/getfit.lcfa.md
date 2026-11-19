# Fit Indices for CFA Models

Fit Indices for CFA Models

## Usage

``` r
# S3 method for class 'lcfa'
getfit(x, digits = 3L, fit_matrix = NULL, ..., model = x)
```

## Arguments

- x:

  Fitted model or collection of fitted models passed to the method.

- digits:

  Number of decimal places used in the returned values.

- fit_matrix:

  Optional precomputed internal fit matrix. This argument is used by
  package post-processing methods to avoid refitting the direct-FIML
  saturated model more than once.

- ...:

  Additional arguments to the generic.

- model:

  A fitted object inheriting from class `"lcfa"`.

## Value

A named numeric vector containing model dimensions and available fit
statistics.
