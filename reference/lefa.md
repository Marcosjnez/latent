# Fit an Exploratory Factor Analysis (EFA) model.

Fit an Exploratory Factor Analysis (EFA) model.

## Usage

``` r
lefa(data, nfactors = 1L, model = NULL, estimator = "ml",
ordered = FALSE, group = NULL,
sample.cov = NULL, nobs = NULL,
positive = FALSE, penalties = TRUE,
missing = "pairwise.complete.obs",
std.lv = FALSE, do.fit = TRUE, mimic = 'latent',
control = NULL, ...)
```

## Arguments

- data:

  data frame or matrix.

- nfactors:

  integer. Number of latent variables.

- estimator:

  Available estimators: "ml", "uls", and "dwls". Defaults to "ml".

- model:

  lavaan's model syntax.

- ordered:

  Logical. Defaults to TRUE.

- group:

  String. Name of the variable that splits the data in different groups.

- sample.cov:

  Covariance matrix between the items. Defaults to NULL.

- nobs:

  Number of observations. Defaults to NULL.

- positive:

  Force a positive-definite solution. Defaults to FALSE.

- penalties:

  list of penalty terms for the parameters.

- missing:

  Method to handle missing data.

- std.lv:

  Provide the parameters of the standardized model.

- do.fit:

  TRUE to fit the model and FALSE to return only the model setup.
  Defaults to TRUE.

- mimic:

  String. Choose the output you want to obtain. Defaults to 'latent'.

- control:

  List of control parameters for the optimization algorithm. See
  'details' for more information.

- ...:

  Additional lavaan arguments. See ?lavaan for more information.

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

`lefa` estimates confirmatory factor models.

## Examples

``` r
if (FALSE) { # \dontrun{
# The famous Holzinger and Swineford (1939) example

fit <- lefa(data = HolzingerSwineford1939, nfactors = 3L)
summary(fit, digits = 3L)
} # }
```
