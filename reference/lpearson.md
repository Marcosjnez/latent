# Pearson Covariance or Correlation Matrix

Estimate a covariance or correlation matrix and its asymptotic
covariance matrix from continuous data.

## Usage

``` r
lpearson(data, model = NULL, std.ov = FALSE,
         VCOV = "standard", likelihood = "normal",
         missing = "pairwise.complete.obs", se = TRUE,
         do.fit = TRUE,
         message = FALSE, control = NULL, ...)
```

## Arguments

- data:

  A data frame or matrix containing numeric observed variables.

- model:

  Optional model specification reserved for internal model setup.

- std.ov:

  Logical. If `TRUE`, the observed variables are standardized before the
  sample matrix is returned, so the off-diagonal elements represent
  correlations rather than covariances.

- VCOV:

  Character string selecting the variance-covariance estimator.
  Available options are `"standard"` and `"robust"`.

- likelihood:

  Character string controlling the covariance denominator. Use
  `"normal"` for the maximum-likelihood denominator \\N\\ and
  `"wishart"` for the usual sample-covariance denominator \\N-1\\.

- missing:

  Character string passed to
  [`stats::cov()`](https://rdrr.io/r/stats/cor.html) to control the
  handling of missing values. `"fiml"` is treated as
  `"pairwise.complete.obs"` because FIML is handled at the CFA level.

- se:

  Logical. Compute and store the covariance matrix of the sample
  covariance/correlation estimates. With `FALSE`, the object is retained
  only as a fixed computational statistic.

- do.fit:

  Logical. If `FALSE`, return the prepared but unfitted `"latent"`
  object.

- message:

  Logical. Print progress messages during estimation.

- control:

  Optional list of internal controls. Custom starting values can be
  supplied through `control$start`.

- ...:

  Additional arguments reserved for future extensions.

## Value

An S4 object of class `"latent"`. The object contains the processed data
in `dataList`, the parameter and optimization structures in `modelInfo`,
the direct estimation output and variance-covariance information in
`Optim`, and the estimated parameter structures in `parameters` and
`transformed_pars`.

## Details

`lpearson()` computes the empirical covariance matrix directly rather
than estimating it numerically. With `std.ov = TRUE`, the sample matrix
is standardized before being stored in the fitted object.

Standard asymptotic covariance matrices are computed under multivariate
normality. With `VCOV = "robust"`, fourth-moment information from the
observed data is used instead.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lpearson(data = HolzingerSwineford1939[, paste0("x", 1:9)])
fit_cor <- lpearson(data = HolzingerSwineford1939[, paste0("x", 1:9)],
                    std.ov = TRUE)
} # }
```
