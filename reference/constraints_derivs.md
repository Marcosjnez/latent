# Derivatives of Parameter Constraints

`constraints_derivs()` is a generic function for obtaining derivatives
of parameter constraints introduced by model transformations.

## Usage

``` r
constraints_derivs(x, ...)
```

## Arguments

- x:

  A fitted model object.

- ...:

  Additional arguments passed to methods.

## Value

A numeric matrix containing derivatives of transformation-induced
parameter constraints.

## Details

Compute derivatives associated with constraints on transformed
parameters.

## See also

[`jacobian`](https://marcosjnez.github.io/latent/reference/jacobian.md),
[`hessian`](https://marcosjnez.github.io/latent/reference/hessian.md),
[`vcov`](https://rdrr.io/r/stats/vcov.html)
