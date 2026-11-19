# Confidence Intervals for Latent Class Models

Compute normal-approximation confidence intervals for parameters of a
fitted latent class model.

## Usage

``` r
# S3 method for class 'llca'
ci(
  x,
  type = "information",
  confidence = 0.95,
  parameters = NULL,
  digits = 3L,
  ...,
  fit = x
)

# S3 method for class 'llca_sam'
ci(x, ..., object = x)

# S3 method for class 'llcalist'
ci(
  x,
  type = "information",
  confidence = 0.95,
  parameters = NULL,
  digits = 3L,
  ...,
  model = x
)
```

## Arguments

- x:

  Fitted model or collection of fitted models passed to the method.

- type:

  Character string selecting the standard-error estimator. Available
  options are `"information"` and `"robust"`. `"standard"` is retained
  as an alias of `"information"`.

- confidence:

  Numeric scalar strictly between zero and one specifying the confidence
  level.

- parameters:

  Optional parameter specification identifying the parameters or
  transformed parameters for which intervals should be returned.

- digits:

  Non-negative integer indicating the number of decimal places used in
  the formatted confidence-interval table.

- ...:

  Additional arguments passed to
  [`se()`](https://marcosjnez.github.io/latent/reference/se.latent.md).

- fit:

  A fitted object inheriting from class `"llca"`.

- object:

  A legacy structural-after-measurement object containing a fitted
  structural `"llca"` component.

- model:

  For the `"llcalist"` method, a collection of fitted `"llca"` models.

## Value

A list containing formatted and numeric confidence limits and the
standard-error result used to construct them.
