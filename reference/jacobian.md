# Jacobian Matrix

`jacobian()` is a generic function for obtaining the derivatives
associated with parameter transformations in fitted latent variable
models.

## Usage

``` r
jacobian(x, ...)
```

## Arguments

- x:

  A fitted model object.

- ...:

  Additional arguments passed to methods.

## Value

A numeric matrix containing derivatives of transformed model parameters.

## Details

Compute Jacobian matrices for transformed model parameters.

## See also

[`hessian`](https://marcosjnez.github.io/latent/reference/hessian.md),
[`vcov`](https://rdrr.io/r/stats/vcov.html),
[`constraints_derivs`](https://marcosjnez.github.io/latent/reference/constraints_derivs.md)
