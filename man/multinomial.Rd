% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multinomial.R
\name{multinomial}
\alias{multinomial}
\title{Family Object for Multinomial Model}
\usage{
multinomial()
}
\value{
An object of class \code{familyLMMsolver} with the following components:
\item{family}{character string with the family name.}
\item{linkfun}{the link function.}
\item{linkinv}{the inverse of the link function.}
\item{dev.resids}{function giving the deviance for each observation as a function of (y, mu, wt)}
}
\description{
The Multinomial model is not part of the standard family. The implementation
is based on Chapter 6 in Fahrmeir et al. (2013).
}
\references{
Fahrmeir, Ludwig, Thomas Kneib, Stefan Lang, Brian Marx, Regression models.
Springer Berlin Heidelberg, 2013.
}
