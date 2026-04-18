# Maximum likelihood estimation of positive-definite polychoric correlation matrices

`lpoly` estimates a polychoric correlation matrix from ordinal data
using maximum likelihood. The function can use either a one-step or a
two-step estimation approach and optionally enforce
positive-semidefiniteness or positive-definiteness through constrained
estimation or penalties.

## Usage

``` r
lpoly(data,
      method = "two-step",
      model = NULL,
      positive = FALSE,
      penalties = FALSE,
      do.fit = TRUE,
      message = FALSE,
      control = NULL,
      ...)
```

## Arguments

- data:

  A data frame or matrix containing the raw ordinal data.

- method:

  Character string indicating the estimation method. Possible values are
  `"one-step"` and `"two-step"`. Default is `"two-step"`.

- model:

  Optional model object used internally for initialization or custom
  model setup. Default is `NULL`.

- positive:

  Logical. If `TRUE`, the estimated polychoric correlation matrix is
  forced to be positive semidefinite. Default is `FALSE`.

- penalties:

  Logical. If `TRUE`, penalties are added to the objective function to
  encourage a positive-definite solution. Default is `FALSE`.

- do.fit:

  Logical. If `TRUE`, the model is fitted. If `FALSE`, only the model
  setup is returned. Default is `TRUE`.

- message:

  Logical. If `TRUE`, progress messages are printed during estimation.
  Default is `FALSE`.

- control:

  A list of control parameters for the optimization algorithm. This may
  include starting values, convergence tolerances, maximum number of
  iterations, and other optimizer-specific options.

- ...:

  Additional arguments passed to internal optimization and model setup
  routines.

## Value

A list containing the fitted model and related information. Typical
elements include:

- `version`: Version number of latent used when the model was estimated.

- `call`: Matched call used to estimate the model.

- `ModelInfo`: Information about the model specification and data.

- `Optim`: Output of the optimization routine.

- `parameters`: Structure containing the model parameters.

- `transparameters`: Structure containing transformed model parameters.

- `loglik`: Log-likelihood of the fitted model.

- `penalized_loglik`: Penalized log-likelihood of the fitted model.

## Details

`lpoly` estimates positive-definite or positive-semidefinite polychoric
correlation matrices from ordinal data. The function is designed for
situations in which the unrestricted polychoric matrix is not guaranteed
to be admissible, for example because of sampling variability or sparse
response patterns.

Two estimation strategies are available:

- `"two-step"`: thresholds are estimated first and the correlation
  matrix is estimated in a second step.

- `"one-step"`: thresholds and correlations are estimated jointly.

If `positive = TRUE`, the estimated matrix is constrained to be positive
semidefinite. If `penalties = TRUE`, penalty terms are added to the
objective function to encourage positive-definiteness.

If `do.fit = FALSE`, the function returns the model setup without
running the optimizer.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lpoly(data = values)

fit_psd <- lpoly(data = values, positive = TRUE)

setup_only <- lpoly(data = values, do.fit = FALSE)
} # }
```
