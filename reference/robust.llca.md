# LatentGold-Style Robust Variance-Covariance Matrix

Construct an information matrix corresponding to a LatentGold-style
sandwich covariance estimator for a fitted latent class model.

## Usage

``` r
# S3 method for class 'llca'
robust(fit)
```

## Arguments

- fit:

  A fitted object of class `"llca"`.

## Value

A symmetric numeric matrix representing the robust information matrix.

## Details

Let \\H\\ denote the Hessian and \\B\\ the empirical covariance matrix
of score contributions. The robust covariance matrix is \$\$H^{-1} B
H^{-1}.\$\$

Because
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md)
expects an information matrix that is subsequently inverted, this
function returns the equivalent matrix \$\$H B^{-1} H.\$\$
