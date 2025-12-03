# Display the sparseness of the mixed model coefficient matrix

Display the sparseness of the mixed model coefficient matrix

## Usage

``` r
displayMME(object, cholesky = FALSE)
```

## Arguments

- object:

  an object of class LMMsolve.

- cholesky:

  Should the cholesky decomposition of the coefficient matrix be
  plotted?

## Value

A plot of the sparseness of the mixed model coefficient matrix.

## Examples

``` r
## Fit model on oats data
data(oats.data)

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = oats.data)

## Obtain deviance.
displayMME(LMM1)

```
