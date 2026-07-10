# Inspect fitted latent class models

Extract model-implied class probabilities, posterior probabilities,
class-conditional response profiles, fit statistics, response-pattern
summaries, classification-error matrices, and optimization diagnostics
from fitted `"llca"` objects.

## Usage

``` r
# S3 method for class 'llca'
latInspect(fit, what = "profile", digits = 4L)

# S3 method for class 'llcalist'
latInspect(model, what = "profile", digits = 4L)
```

## Arguments

- fit:

  An object of class `"llca"` returned by
  [`lca`](https://marcosjnez.github.io/latent/reference/lca.md).

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

  Matching is case-insensitive.

- digits:

  Non-negative integer retained for compatibility with other inspection
  methods. Numeric results are returned without rounding so they can
  safely be used in subsequent computations.

- model:

  An object of class `"llcalist"` containing fitted `"llca"` objects.

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

Gaussian item matrices are read by row name, using `"mean"` and
`"stdv"`. Multinomial item matrices are assumed to contain probabilities
in their first rows and log-probability parameters in their following
rows; only the probability rows are returned in response profiles.

## References

None yet.
