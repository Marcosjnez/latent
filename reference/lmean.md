# Sample Means

Estimate the sample means of a set of observed variables and their
asymptotic covariance matrix.

## Usage

``` r
lmean(data, model = NULL, std.ov = FALSE,
      do.fit = TRUE, message = FALSE,
      control = NULL, ...)
```

## Arguments

- data:

  A data frame or matrix containing numeric observed variables. Column
  names are required and are used as the observed-variable labels.

- model:

  Optional model specification reserved for internal model setup.

- std.ov:

  Logical. If `TRUE`, the observed variables are treated as standardized
  and all mean parameters are fixed to zero. Their asymptotic covariance
  matrix and standard errors are consequently also zero.

- do.fit:

  Logical. If `TRUE`, compute the sample means. If `FALSE`, return the
  prepared but unfitted `"latent"` object.

- message:

  Logical. Print progress messages during estimation.

- control:

  Optional list of internal control parameters. Custom starting values
  can be supplied through `control$start`.

- ...:

  Additional arguments reserved for future extensions.

## Value

An S4 object of class `"latent"`. The object contains the processed data
in `dataList`, the parameter and model structures in `modelInfo`, the
estimated means and standard-error information in `Optim`, and the
parameter-shaped results in `parameters` and `transformed_pars`.

## Details

`lmean()` computes the arithmetic mean of each observed variable using
all available non-missing observations for that variable. When
`std.ov = FALSE`, the variance-covariance matrix of the sample means
retains covariances between variables. With complete data, it is the
sample covariance matrix divided by the sample size. With missing data,
each entry is adjusted for the numbers of observations contributing to
the two means. Standard errors are the square roots of its diagonal.

When `std.ov = TRUE`, the mean parameters are fixed to zero rather than
estimated. This is the appropriate mean structure for standardized
observed variables. Because these values are fixed, their
variance-covariance matrix and standard errors are set to zero.

The function uses the same parameter/model infrastructure as the other
estimators in latent, even though the sample means themselves are
obtained directly rather than through numerical optimization.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lmean(data = HolzingerSwineford1939[, paste0("x", 1:9)])

fit_std <- lmean(data = HolzingerSwineford1939[, paste0("x", 1:9)],
                 std.ov = TRUE)
} # }
```
