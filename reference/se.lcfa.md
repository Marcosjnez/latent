# Standard Errors for Confirmatory Factor Analysis

Extract or compute standard errors for a fitted `"lcfa"` object.

## Usage

``` r
# S3 method for class 'lcfa'
se(fit, type = "standard", digits = 5L, ...)
```

## Arguments

- fit:

  A fitted object of class `"lcfa"`.

- type:

  Character string identifying the requested standard-error type.

- digits:

  Optional number of decimal places used for formatted output.

- ...:

  Additional arguments.

## Value

A list containing the standard errors, variance-covariance matrix,
Hessian, sandwich middle matrix, and the parameter-shaped standard-error
table.
