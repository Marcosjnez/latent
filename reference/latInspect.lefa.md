# Inspect Fitted Exploratory Factor Analysis Objects

Inspect Fitted Exploratory Factor Analysis Objects

## Usage

``` r
# S3 method for class 'lefa'
latInspect(fit, what = "est")
```

## Arguments

- fit:

  A fitted object inheriting from class `"lefa"`.

- what:

  Character string identifying the requested component.

## Value

A rotated parameter list, the unrotated `lcfa` object, or a component
delegated to
[`latInspect.lcfa()`](https://marcosjnez.github.io/latent/reference/latInspect.lcfa.md).
