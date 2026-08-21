# Jacobian Matrix for Latent Models

Compute the cumulative Jacobian from the freely estimated parameters to
selected transformed parameters.

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
  parameters whose derivatives should be returned.

## Value

A numeric matrix whose rows correspond to selected transformed
parameters and whose columns correspond to the freely estimated
parameters.

## Details

The cumulative Jacobian is computed only when this method is called.
Local transformation Jacobians are composed in dependency order, and the
matrix is stored relative to the freely estimated parameter vector
rather than the complete transformed-parameter vector.
