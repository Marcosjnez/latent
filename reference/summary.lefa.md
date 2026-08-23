# Summarize a Fitted Exploratory Factor Analysis

Summarize a Fitted Exploratory Factor Analysis

## Usage

``` r
# S3 method for class 'lefa'
summary(object, digits = 3L, fit = NULL, ...)
```

## Arguments

- object:

  A fitted object inheriting from class `"lefa"`.

- digits:

  Number of decimal places used when printing results.

- fit:

  Alias retained for compatibility. If supplied, `object` is ignored.

- ...:

  Additional arguments reserved for future summary options.

## Value

Invisibly, an object of class `"summary.lefa"` containing the
convergence information, fit indices, and rotated parameter table.
