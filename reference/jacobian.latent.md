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

Only transformations required by the selected outputs are evaluated.
Local Jacobians are composed in dependency order, so chains of
transformations are represented correctly.
