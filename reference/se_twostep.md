# Two-step standard-error adjustment

Adjusts the covariance matrix of a structural model for uncertainty in
the measurement-model parameters estimated in a previous step.

## Usage

``` r
se_twostep(fit2, type = "standard", parameters = NULL)
```

## Arguments

- fit2:

  A fitted structural `"llca"` object whose optimizer control stores the
  fitted measurement model.

- type:

  Character string indicating whether standard or robust covariance
  matrices are used in the two steps.

## Value

A list containing the combined covariance matrix, standard errors, and
the correction matrix `B`.
