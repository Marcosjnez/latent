# Print local bivariate residual diagnostics

Prints the root mean square residual and the pairwise residual summary
table stored in an object returned by `lbvr()`.

## Usage

``` r
# S3 method for class 'lbvr'
print(x, ...)
```

## Arguments

- x:

  An object of class `"lbvr"`.

- ...:

  Additional arguments passed to
  [`print.data.frame()`](https://rdrr.io/r/base/print.dataframe.html)
  when printing the residual summary table.

## Value

The input object `x`, returned invisibly.
