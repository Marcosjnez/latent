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

A list containing the Hessian, variance-covariance matrix, standard
errors, and covariance-method metadata.

## Details

The fitted object is never modified. For ordinary CFA likelihood models,
a temporary derivative copy replaces `cfa_fml` by `cfa_ml` and
`cfa_means_fml` by `cfa_means_ml`. Stored covariance matrices are reused
only when they were explicitly produced by an information/standard
method.
