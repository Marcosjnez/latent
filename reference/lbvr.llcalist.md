# Local bivariate residuals for a collection of latent class models

Applies `lbvr()` to every fitted `"llca"` object contained in an object
of class `"llcalist"`. Model names are preserved. Unnamed models are
labelled according to their number of latent classes.

## Usage

``` r
# S3 method for class 'llcalist'
lbvr(x, digits = 4L, ...)
```

## Arguments

- x:

  An object of class `"llcalist"` containing fitted `"llca"` models.

- digits:

  A non-negative integer indicating the number of decimal places used to
  round the returned diagnostics. The default is `4L`.

- ...:

  Additional arguments passed to the
  [`lbvr.llca()`](https://marcosjnez.github.io/latent/reference/lbvr.llca.md)
  method.

## Value

A named list containing one object of class `"lbvr"` for each fitted
model. The returned list has class `"lbvr.llcalist"`.

## Examples

``` r
if (FALSE) { # \dontrun{
fits <- lca(data = gss82, nclasses = 2:4,
            multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"))
residuals <- lbvr(fits)
} # }
```
