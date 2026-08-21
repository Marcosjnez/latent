# Summarize a Fitted CFA Model

Summarize a Fitted CFA Model

## Usage

``` r
# S3 method for class 'lcfa'
summary(object, digits = 3L, fit = NULL, ...)
```

## Arguments

- object:

  A fitted object inheriting from class `"lcfa"`.

- digits:

  Number of decimal places used when printing results.

- fit:

  Alias retained for backward compatibility. If supplied, `object` is
  ignored.

- ...:

  Additional arguments reserved for future summary options.

## Value

Invisibly, an object of class `"summary.lcfa"` containing the
convergence information, fit indices, and free-parameter table.
