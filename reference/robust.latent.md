# Robust Variance-Covariance Matrix for General Latent Models

Robust Variance-Covariance Matrix for General Latent Models

## Usage

``` r
# S3 method for class 'latent'
robust(fit)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

## Value

A stored robust covariance result when available. If no class-specific
robust covariance exists, the information covariance is returned with a
warning and an explicit fallback label.
