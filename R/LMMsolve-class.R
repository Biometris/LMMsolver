#' Fitted LMMsolve Object
#'
#' An object of class \code{LMMsolve} returned by the LMMsolve function,
#' representing a fitted linear mixed model. Objects of this class have
#' methods for the generic functions coef, .....
#'
#' @return
#' An object of class \code{LMMsolve} contains the following components.
#' \describe{
#' \item{logL}{The restricted log-likelihood at convergence}
#' \item{sigma2e}{The residual error}
#' \item{tau2e}{The estimated variance components}
#' \item{ED}{The effective dimensions}
#' \item{EDmax}{The maximal effective dimensions}
#' \item{EDnames}{The names of the effective dimensions}
#' \item{a}{The estimated effects from the mixed model equations}
#' \item{yhat}{The fitted values}
#' \item{dim}{The dimensions for each of the fixed and random terms in the
#' mixed model}
#' \item{term.labels}{The Names of the fixed and random terms in the mixed
#' model}
#' \item{splRes}{An object with definition of spline argument}
#' }
#'
#' @usage NULL
#'
#' @rdname LMMsolveObject
#' @export
LMMsolveObject <- function(object) {
  structure(object,
            class = c("LMMsolve", "list"))
}
