# Hidden Markov Models

Extend [`lca()`](https://marcosjnez.github.io/latent/reference/lca.md)
with a common transition matrix and covariate-dependent initial state
probabilities. Indicator and outcome distributions, parameter
constraints, starts, penalties, and numerical optimizers use the LCA
helpers.

## Usage

``` r
lhmm(
  data,
  nclasses = 1L,
  gaussian = NULL,
  multinomial = NULL,
  covariates = NULL,
  outcomes = NULL,
  model = NULL,
  weights = NULL,
  adjustment = "bk",
  classification = "modal",
  penalties = TRUE,
  start = NULL,
  control = NULL,
  do.fit = TRUE,
  verbose = TRUE,
  time = NULL,
  id = NULL
)
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
  associations. If an object of class `"llca"` is supplied, its model
  parameters are reused while the class-membership regression
  coefficients are re-estimated. If the supplied object also inherits
  from `"multistage"`, previous `"llca"` objects stored in its `extra`
  slot are traversed recursively so that parameters fixed in earlier
  stages are retained.

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

  Optional list of optimization and estimation controls. Common entries
  include `rstarts` for the number of random starts, `cores` for
  parallel computation, `maxit` for the maximum number of optimizer
  iterations, and `opt` for the optimization algorithm. Available
  optimization methods include `"grad"`, `"lbfgs"`, `"newton"`, and
  `"em"`. Convergence controls include `eps`, `df_eps`, and `step_eps`.
  When `opt = "em"`, `maxit` controls the maximum number of outer EM
  iterations. The numerical optimizer used in the M step is selected
  with `mopt`; available options are `"grad"`, `"lbfgs"`, and
  `"newton"`. `mstep_maxit` controls the maximum number of iterations
  within each M step and `mstep_eps` its convergence tolerance. Missing
  entries are replaced by internal defaults.

- do.fit:

  Logical. If `TRUE`, the model is estimated. If `FALSE`, the function
  returns an `"llca"` object containing the processed data, model
  structure, and optimization setup, but without running the optimizer.

- verbose:

  Logical. If `TRUE`, progress information is printed, especially when
  fitting several values of `nclasses`.

- time:

  Name of the time variable. Observations are sorted within each
  trajectory. With `NULL`, call
  [`lca()`](https://marcosjnez.github.io/latent/reference/lca.md) with
  the same arguments.

- id:

  Name of the trajectory variable. With `NULL`, all observations form
  one trajectory.

## Value

An object of class `"llhmm"`, inheriting from `"latent"`, or an
`"llhmmlist"` for several state counts. Non-temporal calls return the
existing LCA objects. Temporal post-estimation methods are separate.

## Details

In temporal models, covariates at the first observation of each
trajectory predict its initial state probabilities. Later covariate
values are unused. Transition probabilities are common across
trajectories and time; columns correspond to the previous state and rows
to the next state.

Weights contain one positive value per trajectory, in order of first
appearance in the input data. Outcomes contribute at their observed time
points, conditional on the current state. An outcome observed once
should appear only once, with missing values at the other times.

Each adjacent pair of rows represents one transition. Time orders rows;
elapsed-time gaps are not modeled. Repeated response patterns are
retained.

Temporal fitting currently uses one-step estimation. When omitted,
`adjustment` is `"none"` for temporal data and retains the LCA default
otherwise. Temporal EM and multistep adjustments are not implemented. If
every trajectory has one observation, fitting delegates to
[`lca()`](https://marcosjnez.github.io/latent/reference/lca.md).

Use `start$transition` for a positive column-stochastic starting matrix,
or `start$logtransition` for reference logits. The first row of
`logtransition` is fixed to zero by default. Parameter constraints use
`model$logtransition`; other blocks have the same names as in LCA.
