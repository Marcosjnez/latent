# Fitted Latent Class Membership Probabilities

`fitted.llca()` returns the fitted latent class membership probabilities
for the observations retained in an `"llca"` model. These probabilities
depend on the class-membership regression model and its covariates, but
not on the observed indicator responses.

## Usage

``` r
# S3 method for class 'llca'
fitted(object, ...)
```

## Arguments

- object:

  A fitted object of class `"llca"`.

- ...:

  Additional arguments. They are currently ignored.

## Value

A numeric matrix with one row per retained subject and one column per
latent class. Row names correspond to the retained observations and
column names identify the latent classes.

## Details

Extract fitted prior latent class membership probabilities for the
subjects used to estimate a latent class model.

The returned values are the prior class-membership probabilities \\P(C_i
= c \mid x_i)\\. They should not be confused with posterior
probabilities, which additionally condition on the observed indicator
responses and can be obtained with
`latInspect(object, what = "posterior")`.

Class probabilities are stored internally for the unique response and
covariate patterns. This method expands them back to the subject level
using the mapping stored in the fitted object. Consequently, the number
and order of rows correspond to the observations retained in
`object@dataList$data`. Observations removed before estimation because
of missing covariates are not included.

If the model has no covariates, every subject receives the same vector
of estimated class proportions.

## References

None yet.

## See also

[`predict.llca`](https://marcosjnez.github.io/latent/reference/predict.llca.md),
`latInspect`
