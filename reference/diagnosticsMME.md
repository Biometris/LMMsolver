# Give diagnostics for mixed model coefficient matrix C and the cholesky decomposition

Give diagnostics for mixed model coefficient matrix C and the cholesky
decomposition

## Usage

``` r
diagnosticsMME(object)
```

## Arguments

- object:

  an object of class LMMsolve.

## Value

A summary of the mixed model coefficient matrix and its choleski
decomposition.

## Examples

``` r
## Fit model on oats data
data(oats.data)

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = oats.data)

## Obtain deviance.
diagnosticsMME(LMM1)
#> Summary of matrix C 
#> Matrix object of class 'spam' of dimension 26x26,
#>     with 168 (row-wise) nonzero elements.
#>     Density of the matrix is 24.9%.
#> Class 'spam' (32-bit)
#> 
#>  Summary of cholesky decomposition of C 
#> Matrix object of class 'spam' of dimension 26x26,
#>     with 98 (row-wise) nonzero elements.
#>     Density of the matrix is 14.5%.
#> Class 'spam' (32-bit)
```
