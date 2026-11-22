# Inspect fitted latent class and hidden Markov models

Extract model-implied class probabilities, posterior probabilities,
class-conditional response profiles, fit statistics, response-pattern
summaries, classification-error matrices, and optimization diagnostics
from fitted `"llca"` objects.

## Usage

``` r
# S3 method for class 'llca'
latInspect(x, what = "profile", digits = 4L, sort = TRUE, ..., fit = x)

# S3 method for class 'llcalist'
latInspect(x, what = "profile", digits = 4L, sort = TRUE, ..., model = x)
```

## Arguments

- x:

  Fitted model or collection of fitted models passed to the method.

- what:

  Character string indicating the object to extract. Available values
  and aliases include:

  `"profile"`

  :   Posterior class sizes and class-conditional indicator profiles. If
      distal outcomes are present, their class-conditional profiles are
      included in an `outcomes` element.

  `"classconditional"`, `"item"`, `"items"`

  :   Class-conditional indicator parameters. Gaussian indicators
      contain means and standard deviations; multinomial indicators
      contain probabilities only.

  `"outcome"`, `"outcomes"`

  :   Class-conditional distal outcome parameters.

  `"class"`, `"classes"`

  :   Average model-implied prior class probabilities, weighted by
      response-pattern weights.

  `"fullclasses"`

  :   Pattern-specific prior class probabilities.

  `"beta"`, `"coef"`, `"coefs"`

  :   Class-membership regression coefficients.

  `"posterior"`

  :   Case-specific posterior class probabilities.

  `"state"`, `"states"`

  :   Modal posterior class assignments.

  `"respconditional"`

  :   For fully multinomial measurement models, probabilities of
      latent-class membership conditional on each response category.

  `"probcat"`

  :   For fully multinomial measurement models, marginal
      response-category probabilities.

  `"data"`

  :   Processed data stored in the fitted object.

  `"measurement"`

  :   Processed indicator data.

  `"pattern"`, `"patterns"`

  :   Unique model-variable patterns and their observed weights.

  `"table"`, `"summary"`

  :   Response-pattern summary table.

  `"loglik_case"`

  :   Case-specific log-likelihood contributions.

  `"loglik_pattern"`

  :   Pattern-specific weighted log-likelihood contributions.

  `"loss"`, `"loglik"`, `"fit_matrix"`

  :   Optimizer loss, log-likelihood, or the complete matrix of fit
      components.

  `"classification"`

  :   Modal and proportional classification-error matrices. Rows
      represent latent classes and columns represent assigned classes.

  `"convergence"`

  :   Optimizer convergence information.

  `"gradient"`

  :   Gradient diagnostics.

  `"timing"`, `"elapsed"`

  :   Elapsed optimization time.

  `"initial"`, `"transition"`

  :   For a hidden Markov model, the initial probabilities by trajectory
      or the common transition matrix. In the transition matrix, columns
      are previous states and rows are next states.

  `"forward"`, `"filtered"`

  :   Filtered state probabilities conditional on the observations
      through the current time.

  `"backward"`

  :   Scaled likelihood of the future observations for each current
      state. Its value is one at a trajectory's final time. Multiplying
      it by the forward probabilities gives the smoothed probabilities.

  `"smoothed"`, `"posterior"`

  :   For a hidden Markov model, the state probabilities given all
      observations in the trajectory.

  `"posterior_transition"`

  :   Joint probabilities of adjacent states given the whole trajectory,
      returned as an array indexed by next state, previous state, and
      destination observation. Trajectory starts have `NA` slices.

  `"viterbi"`

  :   Most probable joint state sequence for each trajectory. Values are
      integer state identifiers in the fitted model.

  `"loglik_trajectory"`

  :   Unweighted observed-data log-likelihood per trajectory. In a
      temporal model, `"loglik_case"` gives the same result. The sum
      weighted by trajectory weights gives the model log-likelihood.

  Matching is case-insensitive.

- digits:

  Non-negative integer retained for compatibility with other inspection
  methods. Numeric results are returned without rounding so they can
  safely be used in subsequent computations.

- sort:

  Logical. Sort only inspection outputs by decreasing weighted posterior
  size. Defaults to TRUE. Class identities and the original reference
  class are retained, including in modal assignments.

- ...:

  Additional arguments to the generic.

- fit:

  An object of class `"llca"` returned by
  [`lca`](https://marcosjnez.github.io/latent/reference/lca.md).

- model:

  For the `"llcalist"` method, a collection containing `"llca"` and/or
  `"llca_sam"` model fits.

## Value

The object requested through `what`. See the description of `what` for
the possible return types.

## Details

Posterior probabilities and log-likelihood contributions are stored
internally by unique model-variable pattern and expanded back to the
retained observations when case-level results are requested.

Pattern weights are used when calculating average class probabilities,
posterior class sizes, and classification-error matrices. Consequently,
frequency-weighted models produce weighted inspection results.

For hidden Markov models, observation-level results follow the original
input row order. Trajectory-level results follow the order in which IDs
first appeared. `"fullclasses"` gives the initial state probabilities,
and `"classes"` averages them using trajectory weights. The class sizes
in `"profile"` and the default state ordering use smoothed occupancy
across all observations, with each trajectory weight applied to its
observations. `"table"` returns the observations, marginal modal states,
and smoothed probabilities. The modal state at each time can differ from
the Viterbi path. Use `sort = FALSE` to keep the state order from
estimation.

Gaussian item matrices are read by row name, using `"mean"` and
`"stdv"`. Multinomial item matrices are assumed to contain probabilities
in their first rows and log-probability parameters in their following
rows; only the probability rows are returned in response profiles.

## References

None yet.
