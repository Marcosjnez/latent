# Confidence intervals for latent class models

Computes confidence intervals for the parameters of a fitted latent
class model.

## Usage

``` r
# S3 method for class 'llca'
ci(fit, type = "standard", confidence = 0.95, digits = 3L, ...)
```

## Arguments

- fit:

  A fitted object of class `"llca"`.

- type:

  Character string indicating the standard-error estimator used to
  construct the intervals. Available options are `"standard"` and
  `"robust"`. See
  [`se.llca()`](https://marcosjnez.github.io/latent/reference/se.llca.md).

- confidence:

  Numeric scalar strictly between zero and one specifying the confidence
  level.

- digits:

  Non-negative integer indicating the number of decimal places used in
  the formatted confidence-interval table.

- ...:

  Additional arguments passed to other methods.

## Value

A list with the following components:

- `table`:

  A list of formatted parameter tables showing the estimate and its
  confidence interval.

- `lower_table`:

  The lower confidence limits arranged in the model parameter structure.

- `upper_table`:

  The upper confidence limits arranged in the model parameter structure.

- `lower`:

  A named numeric vector of lower confidence limits.

- `upper`:

  A named numeric vector of upper confidence limits.

- `se`:

  A named numeric vector of standard errors.

- `vcov`:

  The variance-covariance matrix used to construct the intervals.

- `B`:

  The empirical or two-step correction matrix, when available.

- `H`:

  The Hessian matrix, when available.

- `newH`:

  The adjusted Hessian used by the robust estimator, when available.

## Details

Confidence limits are computed using the asymptotic normal
approximation. The critical value is obtained as the square root of the
corresponding one-degree- of-freedom chi-squared quantile. Standard
errors are obtained through `se()`, so two-step and robust adjustments
are used when requested.

The `digits` argument affects only the formatted `table`. The numeric
confidence limits, standard errors, and covariance matrices are returned
without rounding.

## See also

`se()`

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lca(data = empathy, nclasses = 3L,
           gaussian = c("ec1", "ec2", "ec3"))

ci(fit)
ci(fit, type = "robust", confidence = 0.90, digits = 4L)
} # }
```
