# Inspect Fitted CFA Objects

Inspect Fitted CFA Objects

## Usage

``` r
# S3 method for class 'lcfa'
latInspect(fit, what = "est")
```

## Arguments

- fit:

  A fitted object inheriting from class `"lcfa"`.

- what:

  Character string identifying the requested component.

## Value

A parameter list, residual list, fit matrix, or estimator-specific
control component, depending on `what`.
