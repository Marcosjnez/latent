# Information Covariance for Multistep Models

Return the multistep covariance matrix implied by the sample-statistic
covariance estimators selected when the model was fitted.

## Usage

``` r
# S3 method for class 'lefa'
information(fit)

# S3 method for class 'multistep'
information(fit)

# S3 method for class 'multistep_lcfa'
information(fit)
```

## Arguments

- fit:

  A fitted object inheriting from class `"multistep"`.

## Value

A list containing the top-level Hessian, covariance matrix, and standard
errors of the freely estimated parameters.
