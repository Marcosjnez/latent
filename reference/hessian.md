# Hessian Matrix

`hessian()` is a generic function for obtaining the matrix of
second-order derivatives of the objective function with respect to the
freely estimated model parameters.

## Usage

``` r
hessian(x, ...)
```

## Arguments

- x:

  A fitted model object.

- ...:

  Additional arguments passed to methods.

## Value

A numeric matrix containing the Hessian of the objective function with
respect to the freely estimated parameters.

## Details

Compute the Hessian matrix of a fitted latent variable model.

## See also

[`vcov`](https://rdrr.io/r/stats/vcov.html),
[`jacobian`](https://marcosjnez.github.io/latent/reference/jacobian.md),
[`constraints_derivs`](https://marcosjnez.github.io/latent/reference/constraints_derivs.md)
