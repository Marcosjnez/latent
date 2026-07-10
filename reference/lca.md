# Latent Class Analysis

`lca()` estimates latent class models in which the measurement model may
contain Gaussian indicators, multinomial indicators, or both. Class
membership probabilities can be unconditional or modeled as a softmax
transformation of observed covariates. Distal outcomes can also be
included and modeled conditionally on the latent classes. The function
supports one-step estimation, Bakk–Kuha two-step estimation, ML
three-step estimation with modal or proportional classification, and
class enumeration by passing several values to `nclasses`.

## Usage

``` r
lca(data, nclasses = 1L, gaussian = NULL, multinomial = NULL,
    covariates = NULL, outcomes = NULL, penalties = TRUE,
    model = NULL, weights = NULL, start = NULL, adjustment = "bk",
    classification = "modal", control = NULL, do.fit = TRUE,
    verbose = TRUE)
```

## Arguments

- data:

  A `data.frame` containing the observed variables. Variables listed in
  `gaussian` and/or `multinomial` are used as indicators. Variables
  listed in `covariates` are used to predict class membership. Variables
  listed in `outcomes` are treated as distal outcomes.

- nclasses:

  A positive integer giving the number of latent classes. If a vector of
  positive integers is supplied, one model is fitted for each value and
  the result is returned as an object of class `"llcalist"`. This is
  usually termed as "class enumeration".

- gaussian:

  Optional character vector with the names of variables in `data` to be
  modeled as Gaussian indicators. These variables are modeled with
  class-specific means and variances. If residual covariance terms
  involving Gaussian indicators are specified in `model`, the
  corresponding variables are modeled jointly with a class-specific
  multivariate Gaussian distribution.

- multinomial:

  Optional character vector with the names of variables in `data` to be
  modeled as multinomial indicators. These variables are converted to
  factors internally. The first factor level is used as the reference
  category for the softmax transformation. Variables with more than 30
  observed categories are not allowed because estimated probabilities
  would be unreliable.

- covariates:

  Optional character vector with the names of variables in `data` used
  to predict latent class membership. Internally, an intercept is always
  included. Character covariates are converted to factors and factors
  are dummy-coded using
  [`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).
  Observations with missing values in the covariates are removed before
  estimation.

- outcomes:

  Optional distal outcomes. It can be either a character vector with
  variable names in `data`, or a named list specifying the likelihood
  for each outcome. If a character vector is supplied, numeric variables
  are treated as Gaussian outcomes and non-numeric variables as
  multinomial outcomes. If a list is supplied, supported names are
  `gaussian` and `multinomial`. For example,
  `outcomes = c("y1", "y2", "z")` automatically assigns numeric outcomes
  to the Gaussian likelihood and non-numeric outcomes to the multinomial
  likelihood; `outcomes = list(gaussian = c("y1", "y2"))` adds Gaussian
  distal outcomes; and `outcomes = list(multinomial = "z")` adds a
  multinomial distal outcome. It is possible to predict different types
  of outcomes at a time using
  `outcomes = list(gaussian = c("y1", "y2"), multinomial = "z")`

- model:

  Optional model specification. This can be a named list used to fix
  parameters, impose equality constraints, or provide custom parameter
  labels. Names should match internal parameter blocks, for example
  `beta`, Gaussian or multinomial item names, or `Sigma|Class<i>` for
  multivariate Gaussian blocks. Character strings containing
  residual-dependency syntax such as `"y1 ~~ y2"` or `"u1 ~~ u2 ~~ u3"`
  are also used to identify residual covariances or residual
  associations. If an object of class `"llca"` is supplied, its
  measurement parameters are reused while the class-membership
  regression coefficients are re-estimated.

- weights:

  Optional numeric vector of observation weights, with one value per row
  of `data`.

- adjustment:

  Character string selecting the estimation strategy when `covariates`
  are supplied. Use `"none"` for one-step estimation, `"bk"` for the
  Bakk–Kuha two-step method, or `"ml"` for the ML three-step correction
  based on classification error. Defaults to `"bk"`.

- classification:

  Character string used when `adjustment = "ml"`. Use `"modal"` for
  modal assignment or `"prop"` for proportional assignment when
  estimating the classification-error matrix.

- penalties:

  Logical value or named list controlling regularization. If `FALSE`, no
  penalties are used. If `TRUE`, default penalties are used. If a named
  list is supplied, missing penalty blocks are filled with their
  defaults. Valid penalty blocks are `beta`, `class`, `prob`, `var`, and
  `Sigma`. The default is
  `list(beta = list(alpha = 0), class = list(alpha = 1), prob = list(alpha = 1), var = list(alpha = 1), Sigma = list(alpha = 1))`.
  The `beta` block may also contain `lambda` and `power` for ridge-type
  regularization.

- start:

  Optional named list of starting values. Names should correspond to
  parameter blocks in the model. Supplied values replace the
  corresponding default initial values, allowing partial specification
  of starting values.

- control:

  Optional list of optimizer and estimation controls. Common entries
  include `rstarts` for the number of random starts, `cores` for
  parallel computation, `maxit` for the maximum number of optimizer
  iterations, `opt` for the optimizer type, and convergence tolerances
  such as `eps`, `df_eps`, and `step_eps`. Missing entries are replaced
  by internal defaults.

- do.fit:

  Logical. If `TRUE`, the model is estimated. If `FALSE`, the function
  returns an `"llca"` object containing the processed data, model
  structure, and optimization setup, but without running the optimizer.

- verbose:

  Logical. If `TRUE`, progress information is printed, especially when
  fitting several values of `nclasses`.

## Value

If `length(nclasses) == 1` and `adjustment = "none"`, an S4 object of
class `"llca"` with the following slots:

- `version`:

  Version of the latent package used to fit the model.

- `call`:

  The matched function call.

- `timing`:

  Elapsed optimization time. Empty when `do.fit = FALSE`.

- `dataList`:

  Processed data objects, including measurement data, covariate design
  matrix, unique response patterns, pattern weights, and mappings
  between original rows and unique patterns.

- `modelInfo`:

  Internal model structure used by the optimizer, including parameter
  labels, transformed-parameter labels, degrees of freedom, manifolds,
  transformations, estimators, and optimizer controls.

- `Optim`:

  Raw output from the optimizer. Empty when `do.fit = FALSE`.

- `parameters`:

  Estimated model parameters on the estimation scale, organized by
  parameter block. Empty when `do.fit = FALSE`.

- `transformed_pars`:

  Estimated parameters after applying model transformations, including
  class probabilities, item probabilities, variances, standard
  deviations, joint probabilities, and log-likelihood components. Empty
  when `do.fit = FALSE`.

- `extra`:

  Additional information reserved for downstream methods.

If `adjustment = "bk"` or `adjustment = "ml"`, the function returns a
list with two elements:

- `measurement`:

  The measurement-model fit from the first step.

- `structural`:

  The structural-model fit with covariates and/or distal outcomes.

If `nclasses` contains several values, a list of fitted objects is
returned with class `"llcalist"`.

## Details

Estimate latent class models with continuous and categorical indicators,
with optional covariates and distal outcomes.

The measurement model is defined by the variables supplied through
`gaussian` and `multinomial`. Variables not listed there are ignored by
the measurement model unless they are used as covariates or distal
outcomes.

For Gaussian indicators, the model estimates class-specific means and
variances. Variances are parameterized through log-variances and
transformed to the positive scale during optimization.

For multinomial indicators, the model estimates class-specific category
probabilities through a softmax parameterization. Each item is stored in
one matrix: the first rows contain probabilities and the following rows
contain their log parameters. The first log parameter is fixed to zero
for identification.

When `covariates` are supplied and `adjustment = "none"`, class
probabilities are modeled in one step through multinomial-log
coefficients stored in the `beta` parameter block. The first class is
the reference class, so its coefficients are fixed to zero.

When `adjustment = "bk"`, the measurement model is first estimated
without covariates. The structural model is then estimated with the
measurement parameters fixed. When `adjustment = "ml"`, the measurement
model is first estimated, cases are assigned to classes using either
modal or proportional classification, and the structural model is
estimated using a classification-error correction.

The function compresses the data into unique response/covariate patterns
and uses their frequencies as pattern weights. If observation weights
are supplied through `weights`, they are incorporated into the pattern
weights.

Residual dependencies can be requested with covariance-style syntax in
`model`. For Gaussian indicators, terms such as `"y1 ~~ y2"` define
class-specific residual covariance parameters. For multinomial
indicators, dependency syntax defines joint categorical response blocks.
For example, `"u1 ~~ u2"` replaces the locally independent likelihood
contribution of `u1` and `u2` by a joint multinomial contribution for
`u1.u2`. Similarly, a connected dependency such as `"u1 ~~ u2 ~~ u3"`
defines one joint variable `u1.u2.u3`. Pairwise interaction parameters
are used to build the joint log-probabilities, and the marginal item
probabilities are recovered from the joint probabilities by summing over
the relevant joint levels.

## References

None yet.

## Examples

``` r
if (FALSE) { # \dontrun{
# Three-class model for categorical indicators
fit <- lca(
  data = gss82,
  nclasses = 3L,
  multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
  penalties = TRUE
)

fit
summary(fit)
getfit(fit)

latInspect(fit, what = "coefs", digits = 3)
latInspect(fit, what = "profile", digits = 3)
latInspect(fit, what = "posterior", digits = 3)

# Class enumeration
fits <- lca(
  data = gss82,
  nclasses = 1:4,
  multinomial = c("PURPOSE", "ACCURACY", "UNDERSTA", "COOPERAT"),
  penalties = TRUE
)

# Latent class regression with Bakk--Kuha adjustment
fit_bk <- lca(
  data = empathy,
  nclasses = 4L,
  gaussian = c("ec1", "ec2", "ec3", "ec4", "ec5", "ec6"),
  covariates = c("sex", "pt1", "pt2", "pt3", "pt4"),
  outcomes = list(gaussian = c("pt5", "pt6")),
  adjustment = "bk"
)
} # }
```
