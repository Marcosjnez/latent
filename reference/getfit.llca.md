# Fit indices for latent class models

Computes likelihood-based fit indices, classification measures, and
entropy for a fitted latent class model or for a list of fitted latent
class models.

## Usage

``` r
# S3 method for class 'llca'
getfit(model, digits = 4)

# S3 method for class 'llcalist'
getfit(model, digits = 4)

# S3 method for class 'llcalist'
getfit(model, digits = 4L)
```

## Arguments

- model:

  An object of class `"llca"` fitted with
  [`lca()`](https://marcosjnez.github.io/latent/reference/lca.md), or an
  object of class `"llcalist"` containing several fitted `"llca"`
  models.

- digits:

  Integer giving the number of decimal places used to round the output.
  Use `NULL` to return the unrounded values.

## Value

For an object of class `"llca"`, a named numeric vector of class
`"getfit.llca"` containing:

- nclasses:

  Number of latent classes.

- npar:

  Number of freely estimated model parameters.

- nobs:

  Number of observations used in estimation.

- loglik:

  Model log-likelihood.

- penalized_loglik:

  Penalized log-likelihood, when penalization is active.

- L2:

  Likelihood-ratio statistic for a fully multinomial model.

- dof:

  Degrees of freedom associated with `L2`.

- pvalue:

  P-value associated with `L2`.

- AIC:

  Akaike information criterion.

- BIC:

  Bayesian information criterion.

- AIC3:

  AIC using a penalty of three per parameter.

- CAIC:

  Consistent Akaike information criterion.

- KIC:

  Kullback information criterion.

- SABIC:

  Sample-size-adjusted BIC.

- ICL:

  Integrated classification likelihood.

- AICp, BICp, AIC3p, CAICp, KICp, SABICp, ICLp:

  Penalized versions of the corresponding indices, returned only when
  penalization is active.

- R2_entropy:

  Entropy-based pseudo-\\R^2\\.

For an object of class `"llcalist"`, a numeric matrix of class
`"getfit.llcalist"` with one row per model and one column per available
fit index.

## Details

For a single `"llca"` model, `getfit()` returns the number of classes,
model parameters, and observations; the log-likelihood; information
criteria; the integrated classification likelihood; and entropy-based
pseudo-\\R^2\\.

The likelihood-ratio statistic \\L^2\\, its degrees of freedom, and its
p-value are returned only when all modeled variables use a multinomial
likelihood. They are returned as `NA` for models containing Gaussian
variables.

When penalization is active, the output also contains the penalized
log-likelihood and the corresponding penalized information criteria. The
attribute `"penalized"` indicates whether these additional indices are
present.

For an `"llcalist"` object, the indices are collected into a matrix with
one row per fitted model.

## References

Akaike, H. (1974). A New Look at the Statistical Model Identification.
*In: Parzen, E., Tanabe, K., Kitagawa, G. (eds) Selected Papers of
Hirotugu Akaike. Springer Series in Statistics. Springer, New York, NY.*
https://doi.org/10.1007/978-1-4612-1694-0_16

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lca(data = gss82, nclasses = 3L,
           multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"))

getfit(fit)
getfit(fit, digits = NULL)

fits <- lca(data = gss82, nclasses = 1:4,
            multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"))

getfit(fits)
} # }
```
