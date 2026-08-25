# Jacobian Matrix for Latent Models

Compute the dependency Jacobian among transformed parameters.

## Usage

``` r
# S3 method for class 'latent'
jacobian(fit, parameters = NULL)
```

## Arguments

- fit:

  A fitted object inheriting from class `"latent"`.

- parameters:

  Optional parameter specification identifying the transformed
  parameters whose Jacobian submatrix should be returned. If `NULL`, the
  complete transformed-parameter Jacobian is returned.

## Value

A sparse square matrix. Rows and columns correspond to the selected
transformed parameters, or to all transformed parameters when
`parameters = NULL`.

## Details

The complete Jacobian has one row and one column for every transformed
parameter. Its diagonal is the identity, while off-diagonal entries
describe direct and transitive dependencies induced by the sequence of
parameter transformations.

This square dependency Jacobian is intended for inspecting relationships
among transformed parameters. It differs from the conventional
delta-method Jacobian, whose columns contain only freely estimated
parameters.
