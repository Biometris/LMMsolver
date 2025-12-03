# Function to obtain restricted log-likelihood and the first derivatives of the log-likelihood, given values for the penalty parameters

Function to obtain restricted log-likelihood and the first derivatives
of the log-likelihood, given values for the penalty parameters

## Usage

``` r
mLogLik(object, theta)
```

## Arguments

- object:

  an object of class LMMsolve

- theta:

  a matrix with values of precision parameters theta.

## Value

A data.frame with logL and the first derivatives of log-likelihood
