# Maximum likelihood estimation of positive definite polychoric correlation matrices.

Maximum likelihood estimation of positive definite polychoric
correlation matrices.

## Usage

``` r
lpoly(data = NULL,
penalties = TRUE,
do.fit = TRUE,
control = NULL)
```

## Arguments

- data:

  data frame or matrix with the raw data.

- penalties:

  Force a positive-definite solution. Defaults to TRUE.

- do.fit:

  TRUE to fit the model and FALSE to return only the model setup.
  Defaults to TRUE.

- control:

  List of control parameters for the optimization algorithm. See
  'details' for more information.

## Value

List with the following objects:

- version:

  Version number of 'latent' when the model was estimated.

- call:

  Code used to estimate the model.

- ModelInfo:

  Model information.

- Optim:

  Output of the optimizer.

- parameters:

  Structure with all model parameters.

- transparameters:

  Structure with all transformed model parameters.

- loglik:

  Logarithm likelihood of the model.

- penalized_loglik:

  Logarithm likelihood + logarithm priors of the model.

## Details

`lpoly` estimates positive-definite polychoric correlation matrices.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lpoly(data = values)
} # }
```
