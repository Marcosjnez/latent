# Derivatives of Parameter Constraints

`constraints_derivs()` is a generic function for obtaining first and
second derivatives of parameter constraints introduced by manifolds and
model transformations.

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

A list containing the first constraint derivatives and the
block-diagonal matrix of second constraint derivatives.

## Details

Compute first and second derivatives associated with constraints on
model parameters and transformed parameters.

## See also

[`jacobian`](https://marcosjnez.github.io/latent/reference/jacobian.md),
[`hessian`](https://marcosjnez.github.io/latent/reference/hessian.md),
[`vcov`](https://rdrr.io/r/stats/vcov.html)
