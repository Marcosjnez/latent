# Fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.

Fit a Confirmatory Factor Analysis (CFA) model with lavaan syntax.

## Usage

``` r
lcfa(data, model = NULL, estimator = "ml",
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

- model:

  lavaan's model syntax.

- estimator:

  Available estimators: "ml", "uls", and "dwls". Defaults to "ml".

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

`lcfa` estimates confirmatory factor models.

## Examples

``` r
if (FALSE) { # \dontrun{
# The famous Holzinger and Swineford (1939) example
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- lcfa(model = HS.model, data = HolzingerSwineford1939)
summary(fit, digits = 3L)
} # }
```
