# Summary of a Latent Class Model

`summary.llca()` displays optimization information, the estimation
method, the number of parameters and observations, available model-fit
statistics, and the estimated latent class profile.

## Usage

``` r
# S3 method for class 'llca'
summary(fit)

# S3 method for class 'llca_sam'
summary(model)

# S3 method for class 'llcalist'
summary(model)
```

## Arguments

- fit:

  A fitted object of class `"llca"`.

- model:

  For the `"llcalist"` method, a collection of fitted `"llca"` and/or
  `"llca_sam"` models.

## Value

The latent class profile returned by
`latInspect(fit, what = "profile")`, invisibly.

## Details

Print a summary of a fitted latent class model.

The printed output reports whether optimization converged, the
optimization method, the number of freely estimated parameters, and
information about the observed response patterns.

For fully multinomial models, the likelihood-ratio statistic, degrees of
freedom, and corresponding p-value are also displayed.

For an `"llca_sam"` object, the summary of the final structural model is
displayed.

For an `"llcalist"`, [`summary()`](https://rdrr.io/r/base/summary.html)
is applied to each model in the collection.

## See also

[`getfit.llca`](https://marcosjnez.github.io/latent/reference/getfit.llca.md),
[`latInspect.llca`](https://marcosjnez.github.io/latent/reference/latInspect.llca.md),
[`se.llca`](https://marcosjnez.github.io/latent/reference/se.llca.md)
