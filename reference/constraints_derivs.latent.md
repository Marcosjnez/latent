# Constraint Derivatives for Latent Models

Some parameter transformations imply constraints on their outputs. For
example, probabilities produced by a softmax transformation sum to one.
`constraints_derivs.latent()` extracts the derivatives associated with
such constraints for selected transformed parameters.

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

A numeric matrix. Rows correspond to the selected transformed parameters
and columns correspond to transformation-induced constraints.

## Details

Compute derivatives of transformation-induced constraints for selected
parameters of a fitted latent variable model.

Only transformations on which the requested parameters depend are
evaluated. Transformations that do not impose an explicit constraint do
not contribute a constraint column.

## See also

[`constraints_derivs`](https://marcosjnez.github.io/latent/reference/constraints_derivs.md),
[`jacobian.latent`](https://marcosjnez.github.io/latent/reference/jacobian.latent.md),
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md)
