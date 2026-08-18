# Robust Variance-Covariance Matrix for Latent Models

Default robust method for latent objects without a class-specific robust
estimator. The information covariance matrix is returned.

## Usage

``` r
# S3 method for class 'latent'
robust(fit)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

## Value

The result of `information(fit)`.
