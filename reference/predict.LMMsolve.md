# Predict function

Predict function

## Usage

``` r
# S3 method for class 'LMMsolve'
predict(
  object,
  newdata,
  type = c("response", "link"),
  se.fit = FALSE,
  deriv = NULL,
  ...
)
```

## Arguments

- object:

  an object of class LMMsolve.

- newdata:

  A data.frame containing new points for which the smooth trend should
  be computed. Column names should include the names used when fitting
  the spline model.

- type:

  When this has the value "link" the linear predictor fitted values or
  predictions (possibly with associated standard errors) are returned.
  When type = "response" (default) fitted values or predictions on the
  scale of the response are returned (possibly with associated standard
  errors).

- se.fit:

  calculate standard errors, default `FALSE`.

- deriv:

  Character string of variable for which to calculate the first
  derivative; default `NULL`.

- ...:

  other arguments. Not yet implemented.

## Value

A data.frame with predictions for the smooth trend on the specified
grid. The standard errors are saved if \`se.fit=TRUE\`.

## Examples

``` r
## simulate some data
f <- function(x) { 0.3 + 0.4*x + 0.2*sin(20*x) }
set.seed(12)
n <- 150
x <- seq(0, 1, length = n)
sigma2e <- 0.04
y <- f(x) + rnorm(n, sd = sqrt(sigma2e))
dat <- data.frame(x, y)

## fit the model
obj <- LMMsolve(fixed = y ~ 1,
         spline = ~spl1D(x, nseg = 50), data = dat)

## make predictions
newdat <- data.frame(x = seq(0, 1, length = 5))
pred <- predict(obj, newdata = newdat, se.fit = TRUE)
pred
#>      x     ypred         se
#> 1 0.00 0.2051355 0.09713998
#> 2 0.25 0.1842741 0.05188443
#> 3 0.50 0.3714877 0.05202181
#> 4 0.75 0.7024734 0.05188443
#> 5 1.00 0.6767744 0.09713998

## make predictions for derivative of x:
pred2 <- predict(obj, newdata = newdat, se.fit = TRUE, deriv = "x")
pred2
#>      x      ypred       se
#> 1 0.00  3.5771184 3.395408
#> 2 0.25  0.5155952 1.549117
#> 3 0.50 -3.5814981 1.474648
#> 4 0.75 -3.1879787 1.549117
#> 5 1.00 -1.0076817 3.395408
```
