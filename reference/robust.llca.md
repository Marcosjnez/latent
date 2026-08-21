# LatentGold-Style Robust Variance-Covariance Matrix

Construct a sandwich covariance estimator for a fitted latent class
model.

## Usage

``` r
# S3 method for class 'llca'
robust(fit)
```

## Arguments

- fit:

  A fitted object of class `"llca"`.

## Value

A list containing the Hessian, empirical score covariance, robust
variance-covariance matrix, standard errors, and method label.
