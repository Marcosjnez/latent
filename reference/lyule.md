# Yule Correlation Matrix

Estimate a matrix of pairwise Yule association coefficients and their
asymptotic standard errors from categorical data.

## Usage

``` r
lyule(data, model = NULL, do.fit = TRUE, control = NULL, ...)
```

## Arguments

- data:

  A data frame or matrix containing categorical observed variables.
  Factor, integer, logical, and numeric variables are supported by the
  underlying Yule association routine.

- model:

  Optional model specification reserved for internal model setup.

- do.fit:

  Logical. If `FALSE`, return the prepared but unfitted `"latent"`
  object.

- control:

  Optional list of internal controls. Custom starting values can be
  supplied through `control$start`.

- ...:

  Additional arguments reserved for future extensions.

## Value

An S4 object of class `"latent"`. The object contains the processed data
in `dataList`, the parameter and optimization structures in `modelInfo`,
the direct estimation output and asymptotic covariance information in
`Optim`, and the estimated parameter structures in `parameters` and
`transformed_pars`.

## Details

`lyule()` estimates the pairwise association matrix directly using the
Yule copula-based association routine implemented in latent. The
diagonal of the resulting matrix is fixed to one and only the
off-diagonal associations are treated as free parameters.

Standard errors are obtained from the pairwise asymptotic standard
errors returned by the Yule association routine. The corresponding
asymptotic covariance matrix is stored in `Optim$SE$ACOV`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lyule(data = data.frame(
  x1 = factor(sample(1:3, 200, replace = TRUE)),
  x2 = factor(sample(1:4, 200, replace = TRUE)),
  x3 = factor(sample(1:2, 200, replace = TRUE))
))
} # }
```
