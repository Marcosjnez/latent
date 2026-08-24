# Robust Variance-Covariance Matrix for Latent Models

Compute or retrieve the robust variance-covariance matrix of the freely
estimated parameters of a fitted latent-variable model.

## Usage

``` r
# S3 method for class 'latent'
robust(fit)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

## Value

A list containing the Hessian when available, the empirical score
covariance when available, the variance-covariance matrix, standard
errors, and covariance-method metadata.

## Details

Deterministic multistep models use the same propagated covariance for
`information()` and `robust()` because the covariance method of every
preceding statistic is fixed when the multistep model is fitted.

Pearson covariance and correlation estimators compute their robust
covariance directly from the observed data. Other latent models reuse a
stored robust covariance when one is available. If no robust covariance
has been defined, the information covariance is returned with a warning.

Direct robust likelihood scores are not yet available for ordinary
`"lcfa"` likelihood models.
