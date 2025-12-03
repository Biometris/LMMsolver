# Summarize Linear Mixed Model fits

Summary method for class "LMMsolve". Creates either a table of effective
dimensions (which = "dimensions") or a table of variances (which =
"variances").

## Usage

``` r
# S3 method for class 'LMMsolve'
summary(object, which = c("dimensions", "variances"), ...)

# S3 method for class 'summary.LMMsolve'
print(x, ...)
```

## Arguments

- object:

  An object of class LMMsolve

- which:

  A character string indicating which summary table should be created.

- ...:

  Some methods for this generic require additional arguments. None are
  used in this method.

- x:

  An object of class summary.LMMsolve, the result of a call to
  summary.LMM

## Value

A data.frame with either effective dimensions or variances depending on
which.

## Methods (by generic)

- `print(summary.LMMsolve)`: print summary

## Examples

``` r
## Fit model on oats data.
data(oats.data)

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = oats.data)

## Obtain table of effective dimensions.
summ1 <- summary(LMM1)
print(summ1)
#> Table with effective dimensions and penalties: 
#> 
#>         Term Effective Model Nominal Ratio Penalty
#>  (Intercept)         1     1       1     1    0.00
#>          rep         2     2       2     1    0.00
#>          gen        23    23      23     1    0.00
#>     residual        46    72      46     1    7.43
#> 
#>  Total Effective Dimension: 72 

## Obtain table of variances.
summ2 <- summary(LMM1,
                which = "variances")
print(summ2)
#> Table with variances: 
#> 
#>   VarComp Variance
#>  residual     0.13
#> 
```
