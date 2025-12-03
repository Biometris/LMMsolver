# Log-likelihood of an LMMsolve object

Obtain the Restricted Maximum Log-Likelihood of a model fitted using
LMMsolve.

## Usage

``` r
# S3 method for class 'LMMsolve'
logLik(object, includeConstant = TRUE, ...)
```

## Arguments

- object:

  an object of class LMMsolve

- includeConstant:

  Should the constant in the restricted log-likelihood be included.
  Default is `TRUE`, as for example in `lme4` and SAS. In `asreml` the
  constant is omitted.

- ...:

  some methods for this generic require additional arguments. None are
  used in this method.

## Value

The restricted maximum log-likelihood of the fitted model.

## Examples

``` r
## Fit model on oats data
data(oats.data)

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = oats.data)

## Obtain log-likelihood.
logLik(LMM1)
#> [1] -34.95557

## Obtain log-likelihood without constant.
logLik(LMM1, includeConstant = FALSE)
#> [1] 7.315605
```
