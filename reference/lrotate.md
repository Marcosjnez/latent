# Rotate the lambda matrix of an orthogonal factor model.

Rotate the lambda matrix of an orthogonal factor model.

## Usage

``` r
lrotate(lambda, projection = "oblq", rotation = "oblimin",
 group = NULL, positive = FALSE, penalties = TRUE,
 do.fit = TRUE, control = NULL, ...)
```

## Arguments

- lambda:

  List, loading matrices for each group.

- projection:

  String. Can be "orth", "oblq", or "poblq".

- rotation:

  String. Name of the variable that splits the data in different groups.

- group:

  String. Name of the variable that splits the data in different groups.

- positive:

  Force a positive-definite solution. Defaults to FALSE.

- penalties:

  list of penalty terms for the parameters.

- do.fit:

  TRUE to fit the model and FALSE to return only the model setup.
  Defaults to TRUE.

- control:

  List of control parameters for the optimization algorithm. See
  'details' for more information.

- ...:

  Additional arguments.

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

`lrotate` estimates confirmatory factor models.

## Examples

``` r
if (FALSE) { # \dontrun{

fit <- lrotate(lambda = , projection = "oblq", rotation = "oblimin")
summary(fit, digits = 3L)
} # }
```
