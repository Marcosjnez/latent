# Predict latent class membership probabilities

\`predict.llca()\` returns prior latent class membership probabilities
implied by the fitted class-membership model. Predictions depend only on
the covariates and the estimated \`beta\` coefficients. Measurement
indicators and distal outcomes are not used.

\`predict.llcalist()\` applies the prediction method to collections of
latent class models. When the object contains named measurement and
structural fits, only predictions from the structural fit are returned.

## Usage

``` r
# S3 method for class 'llca'
predict(model, new = NULL, ...)

# S3 method for class 'llcalist'
predict(model, new = NULL, ...)

# S3 method for class 'llcalist'
predict(model, new = NULL, ...)
```

## Arguments

- model:

  A fitted object of class \`"llca"\` or \`"llcalist"\`.

- new:

  Optional \`data.frame\` or matrix containing new values for all
  covariates used in the fitted model. Column names must match the
  original covariate names. When \`NULL\`, predictions are computed for
  the data used to fit the model. For a model without covariates,
  \`new\` only determines the number and row names of the returned
  predictions.

- ...:

  Additional arguments. Currently not used.

## Value

For an \`"llca"\` object, a \`data.frame\` containing the non-intercept
columns of the prediction design matrix, followed by one latent class
probability column per class. If the model has no covariates, only the
class probability columns are returned.

For an \`"llcalist"\`, a single \`data.frame\` is returned when only one
relevant structural model is present. Otherwise, a list of prediction
\`data.frame\`s is returned with class \`"predict.llcalist"\`.

## Details

Compute latent class membership probabilities from a fitted latent class
model, using the class-membership regression coefficients and, when
present, observed covariates.

For models with covariates, the method reconstructs the design matrix
using the variable types and factor levels from the fitted data. Its
columns are then matched to the design matrix used during estimation
before calculating the linear predictors and applying the softmax
transformation.

For models without covariates, the fitted latent class probabilities are
unconditional. The same probability vector is returned for every
original observation, or for every row of \`new\` when new data are
supplied.

The returned predictions are prior class-membership probabilities based
only on the class-membership model. They are different from posterior
class probabilities, which also condition on the observed indicators and
can be obtained with \`latInspect(model, what = "posterior")\`.

For an \`"llcalist"\` containing named measurement and structural fits,
the measurement fit is skipped. For other collections containing both
models with and without covariates, only models with covariates are
retained. For class-enumeration objects, predictions are returned for
every relevant fitted model.

## References

None yet.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lca(data = gss82, nclasses = 3L,
           multinomial = c("PURPOSE", "ACCURACY"))

predict(fit)

fit_cov <- lca(data = empathy, nclasses = 3L,
               gaussian = c("ec1", "ec2", "ec3"),
               covariates = c("sex", "pt1"), adjustment = "none")

predict(fit_cov)
predict(fit_cov, new = data.frame(sex = c("female", "male"),
                                  pt1 = c(2.1, 3.4)))
} # }
```
