# Robust Covariance for Multistep Models

The covariance estimators of preceding sample statistics are selected
when the multistep model is fitted. This method therefore returns the
same multistep propagation as
[`information.multistep()`](https://marcosjnez.github.io/latent/reference/information.multistep.md).

## Usage

``` r
# S3 method for class 'lefa'
robust(fit)

# S3 method for class 'multistep'
robust(fit)

# S3 method for class 'multistep_lcfa'
robust(fit)
```

## Arguments

- fit:

  A fitted object inheriting from class `"multistep"`.

## Value

The result of `information.multistep(fit)`.
