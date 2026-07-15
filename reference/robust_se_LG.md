# LatentGold-style robust standard errors

Computes a sandwich covariance estimator from the Hessian and the score
contribution of each observed response pattern.

## Usage

``` r
robust_se_LG(fit, parameters = NULL)
```

## Arguments

- fit:

  A fitted object of class `"llca"`.

## Value

A list containing the standard errors, covariance matrix, empirical
score covariance matrix `B`, Hessian `H`, and adjusted Hessian `newH`.
