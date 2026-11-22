# Inspect Fitted CFA and SEM Objects

Inspect Fitted CFA and SEM Objects

## Usage

``` r
# S3 method for class 'lcfa'
latInspect(x, what = "est", sort = TRUE, ..., fit = x)
```

## Arguments

- x:

  Fitted model or collection of fitted models passed to the method.

- what:

  Character string identifying the requested component. Use `"est"` or
  `"labels"` for the model parameters, or request a block such as
  `"lambda"`, `"theta"`, `"psi"` or `"intercepts"`. Use `"implied"`,
  `"implied.means"`, `"implied.thresholds"` or `"residuals"` for model
  moments. `"loglik"` and `"loss"` return fit components by group and
  overall. Log-likelihoods are available for likelihood estimators. Use
  `"convergence"`, `"gradient"` or `"elapsed"` for optimizer
  diagnostics. Elapsed time is in seconds.

- sort:

  Logical. Sort and orient factor-specific outputs. Defaults to TRUE.

- ...:

  Additional arguments to the generic.

- fit:

  A fitted object inheriting from class `"lcfa"`.

## Value

A parameter list, residual list, fit matrix, or estimator-specific
control component, depending on `what`.
