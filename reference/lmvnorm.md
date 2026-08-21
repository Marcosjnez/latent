# Saturated Multivariate-Normal Moments with Incomplete Data

Estimate an unrestricted vector of means and variance-covariance matrix
by maximizing the multivariate-normal observed-data likelihood.
Missing-data patterns are represented internally by their means,
covariance matrices, and frequencies, but only one global mean vector
and covariance matrix are estimated per substantive group.

## Usage

``` r
lmvnorm(data, group = NULL, variables = NULL,
        se = TRUE, do.fit = TRUE, message = FALSE,
        control = NULL, ...)
```

## Arguments

- data:

  A data frame or numeric matrix containing the observed variables.

- group:

  Optional character string identifying a grouping variable.

- variables:

  Optional character vector of observed-variable names, or a list
  containing one character vector per group. If `NULL`, all numeric
  columns other than the grouping variable are used.

- se:

  Logical or character. `TRUE`, `"standard"`, and `"information"`
  compute the information covariance matrix. `FALSE` skips its
  computation. Robust covariance estimation is not yet implemented for
  this estimator.

- do.fit:

  Logical. If `FALSE`, return the prepared but unfitted `"latent"`
  object.

- message:

  Logical. Print progress messages.

- control:

  Optional list of optimizer controls. Custom starting values may be
  supplied through `control$start`.

- ...:

  Additional arguments reserved for future extensions.

## Value

An S4 object of class `"latent"`. Its public parameter blocks are named
`means` and `S` for a single group, with group suffixes for multiple
groups. The joint finite-sample covariance of all estimated moments is
stored in `Optim$SE$VCOV`.

## Details

For every missingness pattern \\r\\, the likelihood contribution is
evaluated from the pattern frequency, mean vector, and covariance
matrix. The covariance matrix uses the maximum-likelihood divisor
\\n_r\\; for a singleton pattern it is the zero matrix.

The optimized objective is the total negative observed-data
log-likelihood. Pattern means and covariance matrices are fixed
transformed parameters, so the information covariance is the inverse
Hessian on its natural finite-sample scale.

## Examples

``` r
if (FALSE) { # \dontrun{
X <- HolzingerSwineford1939[, paste0("x", 1:9)]
X$x1[1:20] <- NA
fit <- lmvnorm(X)
fit@transformed_pars$means
fit@transformed_pars$S
} # }
```
