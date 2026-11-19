# Summary of a Latent Class Model

Print optimization information, likelihood-based fit indices, and the
estimated latent class profile of a fitted latent class model.

## Usage

``` r
# S3 method for class 'llca'
summary(object, digits = 3L, ..., fit = object)

# S3 method for class 'llcalist'
summary(object, digits = 3L, ..., model = object)
```

## Arguments

- object:

  Fitted model or collection of fitted models passed to the method.

- digits:

  Non-negative integer giving the number of decimal places used in
  printed numeric results.

- ...:

  Additional arguments reserved for future summary options.

- fit:

  A fitted object inheriting from class `"llca"`.

- model:

  For the `"llcalist"` method, a collection of fitted `"llca"` models.

## Value

The latent class profile returned by
`latInspect(fit, what = "profile")`, invisibly.
