# Fitted Probabilities for Lists of Latent Class Models

Extract fitted prior latent class membership probabilities from an
`"llcalist"` object.

## Usage

``` r
# S3 method for class 'llcalist'
fitted(object, ...)
```

## Arguments

- object:

  An object of class `"llcalist"` containing fitted `"llca"` and/or
  `"llca_sam"` objects.

- ...:

  Additional arguments passed to the corresponding
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) methods.

## Value

A named list containing fitted latent class membership probabilities for
each model in `object`. The result has class `"fitted.llcalist"`.

## See also

[`fitted.llca`](https://marcosjnez.github.io/latent/reference/fitted.llca.md),
[`predict.llca`](https://marcosjnez.github.io/latent/reference/predict.llca.md)
