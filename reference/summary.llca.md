# Summary of a Latent Class Model

Print optimization information, likelihood-based fit indices, and the
estimated latent class profile of a fitted latent class model.

## Usage

``` r
# S3 method for class 'llca'
summary(fit, digits = 3L, ...)

# S3 method for class 'llcalist'
summary(model, digits = 3L, ...)
```

## Arguments

- fit:

  A fitted object inheriting from class `"llca"`.

- digits:

  Non-negative integer giving the number of decimal places used in
  printed numeric results.

- ...:

  Additional arguments reserved for future summary options.

- model:

  For the `"llcalist"` method, a collection of fitted `"llca"` models.

## Value

The latent class profile returned by
`latInspect(fit, what = "profile")`, invisibly.
