# Inspect Fitted Exploratory Factor Analysis Objects

Inspect Fitted Exploratory Factor Analysis Objects

## Usage

``` r
# S3 method for class 'lefa'
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

  A fitted object inheriting from class `"lefa"`.

## Value

A rotated parameter list, the unrotated `lcfa` object, or a component
delegated to
[`latInspect.lcfa()`](https://marcosjnez.github.io/latent/reference/latInspect.lcfa.md).
