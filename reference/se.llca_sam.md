# Standard Errors for Structural-After-Measurement Models

Compute standard errors for the structural component of a multi-step
latent class model.

## Usage

``` r
# S3 method for class 'llca_sam'
se(fit, type = "standard", parameters = NULL, digits = 4L, ...)
```

## Arguments

- fit:

  An object of class `"llca_sam"` containing fitted `measurement` and
  `structural` components.

- type:

  Character string specifying the covariance estimator. See
  [`se.llca`](https://marcosjnez.github.io/latent/reference/se.llca.md).

- parameters:

  Optional parameter specification passed to
  [`se.llca`](https://marcosjnez.github.io/latent/reference/se.llca.md).

- digits:

  Non-negative integer specifying the number of decimal places used in
  formatted parameter tables. Use `NULL` to avoid rounding.

- ...:

  Additional arguments passed to
  [`se.llca`](https://marcosjnez.github.io/latent/reference/se.llca.md).

## Value

The result of `se(fit$structural, ...)`.

## Details

The method delegates inference to the structural `"llca"` model. When
that structural model stores the fitted measurement model among its
model specifications,
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md)
automatically propagates measurement-model uncertainty to the structural
estimates.

## See also

[`se.llca`](https://marcosjnez.github.io/latent/reference/se.llca.md),
[`vcov.latent`](https://marcosjnez.github.io/latent/reference/vcov.latent.md),
[`lca`](https://marcosjnez.github.io/latent/reference/lca.md)
