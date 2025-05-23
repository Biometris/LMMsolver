% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/displayMME.R
\name{displayMME}
\alias{displayMME}
\title{Display the sparseness of the mixed model coefficient matrix}
\usage{
displayMME(object, cholesky = FALSE)
}
\arguments{
\item{object}{an object of class LMMsolve.}

\item{cholesky}{Should the cholesky decomposition of the coefficient matrix
be plotted?}
}
\value{
A plot of the sparseness of the mixed model coefficient matrix.
}
\description{
Display the sparseness of the mixed model coefficient matrix
}
\examples{
## Fit model on john.alpha data from agridat package.
data(john.alpha, package = "agridat")

## Fit simple model with only fixed effects.
LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
                data = john.alpha)

## Obtain deviance.
displayMME(LMM1)

}
