# LatentGold-style residual for two continuous indicators

Computes the score-based raw bivariate residual used by `lbvr()` for a
pair of continuous indicators. The fitted model is recreated with the
selected residual dependency while all previously estimated parameters
remain fixed.

## Usage

``` r
latentgold_cont_cont_bvr(fit, v1, v2)
```

## Arguments

- fit:

  A fitted object of class `"llca"`.

- v1:

  Character string naming the first continuous indicator.

- v2:

  Character string naming the second continuous indicator.

## Value

A list with components `resid`, `pval`, `score_stat`, and `df`.
