% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LMMsolve-class.R
\name{residuals.LMMsolve}
\alias{residuals.LMMsolve}
\title{Residuals of an LMMsolve object.}
\usage{
\method{residuals}{LMMsolve}(object, ...)
}
\arguments{
\item{object}{an object of class LMMsolve}

\item{...}{some methods for this generic require additional arguments.
None are used in this method.}
}
\value{
A vector of residuals.
}
\description{
Obtain the residuals from a mixed model fitted using LMMSolve.
}
\examples{
## Fit model on john.alpha data from agridat package.
data(john.alpha, package = "agridat")

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = john.alpha)

## Obtain fitted values.
residuals1 <- residuals(LMM1)

}
