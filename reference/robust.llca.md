# LatentGold-Style Robust Variance-Covariance Matrix

Construct a sandwich covariance estimator for a fitted latent class
model.

## Usage

``` r
# S3 method for class 'llca'
robust(x, ..., fit = x)
```

## Arguments

- x:

  Fitted model or collection of fitted models passed to the method.

- ...:

  Additional arguments to the generic.

- fit:

  A fitted object of class `"llca"`.

## Value

A list containing the Hessian, empirical score covariance, robust
variance-covariance matrix, standard errors, and method label.
