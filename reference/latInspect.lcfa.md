# Inspect Fitted CFA and SEM Objects

Inspect Fitted CFA and SEM Objects

## Usage

``` r
# S3 method for class 'lcfa'
latInspect(x, what = "est", sort = TRUE, ..., fit = x)
```

## Arguments

- x:

  Fitted model or collection of fitted models passed to the method.

- what:

  Character string identifying the requested component.

- sort:

  Logical. Sort/orient only the requested output; defaults to TRUE.

- ...:

  Additional arguments to the generic.

- fit:

  A fitted object inheriting from class `"lcfa"`.

## Value

A parameter list, residual list, fit matrix, or estimator-specific
control component, depending on `what`.
