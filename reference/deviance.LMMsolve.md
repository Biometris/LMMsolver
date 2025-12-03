# Deviance of an LMMsolve object

Obtain the deviance of a model fitted using LMMsolve.

## Usage

``` r
# S3 method for class 'LMMsolve'
deviance(object, relative = TRUE, includeConstant = TRUE, ...)
```

## Arguments

- object:

  an object of class LMMsolve

- relative:

  Deviance relative conditional or absolute unconditional
  (-2\*logLik(object))? Default `relative = TRUE`.

- includeConstant:

  Should the constant in the restricted log-likelihood be included.
  Default is `TRUE`, as for example in `lme4` and SAS. In `asreml` the
  constant is omitted.

- ...:

  some methods for this generic require additional arguments. None are
  used in this method.

## Value

The deviance of the fitted model.

## Examples

``` r
## Fit model on oats.data
data(oats.data)

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = oats.data)

## Obtain deviance.
deviance(LMM1)
#> [1] 6.190954
```
