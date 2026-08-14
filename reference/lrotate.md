# Rotate a factor loading matrix

`lrotate` rotates one or more factor loading matrices using an
orthogonal or oblique projection and a selected rotation criterion.

## Usage

``` r
lrotate(lambda, projection = "oblq", rotation = "oblimin",
        group = NULL, positive = FALSE, penalties = TRUE,
        do.fit = TRUE, control = NULL, ...)
```

## Arguments

- lambda:

  A matrix or a list of loading matrices, one for each group.

- projection:

  Character string. Available projections are `"orth"`, `"oblq"`, and
  `"poblq"`.

- rotation:

  Character string identifying the rotation criterion.

- group:

  Optional grouping information retained in the fitted object.

- positive:

  Logical. Retained for compatibility with the factor-analysis
  interface.

- penalties:

  Logical or list of penalty settings retained in the optimization
  control.

- do.fit:

  Logical. If `TRUE`, fit the rotation. If `FALSE`, return only the
  model specification.

- control:

  List of optimization-control arguments.

- ...:

  Additional arguments required by the selected projection or rotation
  criterion.

## Value

An object of class `"latent"`.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- lrotate(lambda = list(lambda),
               projection = "oblq",
               rotation = "oblimin")
} # }
```
