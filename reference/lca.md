# Latent Class Analysis.

Estimate latent class models with gaussian and multinomial item models,
with or without covariates to predict class membership.

## Usage

``` r
lca(data, nclasses = 2L, item = rep("gaussian", ncol(data)),
    X = NULL, penalties = TRUE, model = NULL, mimic = "LG",
    start = NULL, do.fit = TRUE, verbose = TRUE, control = NULL)
```

## Arguments

- data:

  data frame or matrix.

- nclasses:

  Number of latent classes.

- item:

  Character vector with the model for each item (i.e., "gaussian" or
  "multinomial"). Defaults to "gaussian" for all the items.

- X:

  Matrix of covariates.

- penalties:

  Boolean or list of penalty terms for the parameters.

- model:

  List of parameter labels. See 'details' for more information.

- mimic:

  String. Replicate the output of other softwares. Use "LG" to replicate
  the output of LatentGOLD.

- start:

  List of starting values for the parameters. See 'details' for more
  information.

- do.fit:

  TRUE to fit the model and FALSE to return only the model setup.
  Defaults to TRUE.

- control:

  List of control parameters for the optimization algorithm. See
  'details' for more information.

- verbose:

  Print information of model estimation. Defaults to FALSE.

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

- user_model:

  Structure with most relevant parameters.

- parameters:

  Structure with all model parameters.

- transparameters:

  Structure with all transformed model parameters.

- posterior:

  Posterior probability of class membership.

- state:

  Class with highest posterior probability.

- loglik:

  Logarithm likelihood of the model.

- penalized_loglik:

  Logarithm likelihood + logarithm priors of the model.

- loglik_case:

  Logarithm likelihood of each pattern.

- summary_table:

  Table of summaries of the fitted model. Useful when all the items are
  'multinomial'.

- ClassConditional:

  Parameters that are conditional on the class memberships.

- ResponseConditional:

  Probability of class memberships conditional on item response. Only
  available when all the items have a 'multinomial' likelihood.

- probCat:

  Marginal probability of an item response. Only available when all the
  items have a 'multinomial' likelihood.

## Details

`lca` estimates models with categorical and continuous data.

## References

None yet.

## Examples

``` r
if (FALSE) { # \dontrun{

fit <- lca(data = gss82, nclasses = 3L,
           item = rep("multinomial", ncol(gss82)),
           penalties = TRUE, do.fit = TRUE)
fit@timing
fit@loglik # -2754.643
fit@penalized_loglik # -2759.507
fit@Optim$opt$iterations

# Plot model fit info:
fit

# Get fit indices:
getfit(fit)

# Get a summary:
summary(fit)

# Inspect model objects:
latInspect(fit, what = "coefs", digits = 3)
latInspect(fit, what = "classes", digits = 3)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)
latInspect(fit, what = "table", digits = 3)
latInspect(fit, what = "pattern", digits = 3)

# Get standard errors:
SE <- se(fit, type = "standard", digits = 4)
SE$table

# Get confidence intervals:
CI <- ci(fit, type = "standard", confidence = 0.95, digits = 2)
CI$table
} # }
```
