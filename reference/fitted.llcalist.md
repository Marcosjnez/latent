# Fitted Probabilities for Lists of Latent Class Models

`fitted.llcalist()` applies
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) to latent class
models stored in an `"llcalist"`. For Bakk–Kuha or ML adjustment results
containing elements named `measurement` and `structural`, only the
structural model is returned because it contains the class-membership
regression.

## Usage

``` r
# S3 method for class 'llcalist'
fitted(object, ...)
```

## Arguments

- object:

  An object of class `"llcalist"` containing fitted `"llca"` objects.

- ...:

  Additional arguments passed to
  [`fitted.llca()`](https://marcosjnez.github.io/latent/reference/fitted.llca.md).

## Value

For an adjusted model with a named `structural` component, a numeric
matrix containing the fitted class-membership probabilities from that
model. Otherwise, a list containing one such matrix for each fitted
`"llca"` model. The latter result has class `"fitted.llcalist"`.

## Details

Extract fitted prior class-membership probabilities from an `"llcalist"`
object.

## References

None yet.

## See also

[`fitted.llca`](https://marcosjnez.github.io/latent/reference/fitted.llca.md),
[`predict.llcalist`](https://marcosjnez.github.io/latent/reference/predict.llca.md)
