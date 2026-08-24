# Constraint Derivatives for Latent Models

Some manifolds and parameter transformations imply constraints. For
example, parameters on a unit manifold have unit norm and probabilities
produced by a softmax transformation sum to one.
`constraints_derivs.latent()` extracts the first and second derivatives
associated with such constraints for selected transformed parameters.

## Usage

``` r
# S3 method for class 'latent'
constraints_derivs(fit, parameters = NULL)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

- parameters:

  Optional parameter specification identifying the transformed
  parameters for which constraint derivatives should be returned.
  Parameter labels must occur in `fit@modelInfo$transparameters_labels`.
  If `NULL`, the parameter blocks corresponding to the fitted model
  parameters are used.

## Value

A list with two sparse matrices:

- `dconstr`:

  Rows correspond to selected transformed parameters and columns
  correspond to constraints.

- `d2constr`:

  A block-diagonal matrix containing one square constraint Hessian per
  constraint, in the same order as the columns of `dconstr`.

## Details

Compute first and second derivatives of manifold- and
transformation-induced constraints for selected parameters of a fitted
latent variable model.

Only transformations on which the requested parameters depend are
evaluated. Manifold constraints are always evaluated. Transformations
and manifolds that do not impose an explicit constraint do not
contribute a constraint column.

If \\p\\ parameters and \\m\\ constraints are selected, `d2constr` has
dimension \\pm\\ by \\pm\\. Its \\j\\-th \\p\\ by \\p\\ diagonal block
is the Hessian of constraint \\j\\.

## See also

[`constraints_derivs`](https://marcosjnez.github.io/latent/reference/constraints_derivs.md),
[`jacobian.latent`](https://marcosjnez.github.io/latent/reference/jacobian.latent.md),
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md)
