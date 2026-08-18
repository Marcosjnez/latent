# Variance-Covariance Matrix for Confirmatory Factor Models

Compute the covariance matrix of CFA parameters by propagating the
sampling covariance matrix of the means, covariances, correlations, and
thresholds used as sample statistics.

## Usage

``` r
# S3 method for class 'lcfa'
robust(fit)
```

## Arguments

- fit:

  A fitted object of class `"lcfa"`.

## Value

A list containing the Hessian, cross-derivative matrix, sandwich middle
matrix, parameter variance-covariance matrix, and standard errors.
