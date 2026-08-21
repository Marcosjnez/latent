# Information Variance-Covariance Matrix for Latent Models

Information Variance-Covariance Matrix for Latent Models

## Usage

``` r
# S3 method for class 'latent'
information(fit)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

## Value

A list containing the Hessian, variance-covariance matrix, and standard
errors of the freely estimated parameters.

## Details

The fitted object is never modified. Before evaluating a Hessian for an
ordinary latent object, a temporary derivative copy replaces `cfa_fml`
by `cfa_ml` and `cfa_means_fml` by `cfa_means_ml`. The latter estimators
evaluate the total negative log-likelihood and share the same
parameter-index interfaces as their FML counterparts. Multistep objects
continue to use their fitted discrepancy functions because their
covariance is propagated by
[`se.multistep()`](https://marcosjnez.github.io/latent/reference/se.multistep.md).
