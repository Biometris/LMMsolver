#' Fitted LMMsolve Object
#'
#' An object of class \code{LMMsolve} returned by the LMMsolve function,
#' representing a fitted linear mixed model. Objects of this class have
#' methods for the generic functions coef, .....
#'
#' @return
#' An object of class \code{LMMsolve} contains the following components:
#' \item{logL}{The restricted log-likelihood at convergence}
#' \item{dev}{The REML deviance at convergence (i.e., - 2 times the
#' restricted log-likelihood \code{logL})}
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
#'
#' @usage NULL
#'
#' @rdname LMMsolveObject
#' @export
LMMsolveObject <- function(object) {
  structure(object,
            class = c("LMMsolve", "list"))
}


#' Summarize Linear Mixed Model fits
#'
#' Summary method for class "LMMsolve". Prints the effective dimensions of the
#' model.
#'
#' @param object an object of class LMMsolve
#' @param which default dimensions (only option so far).
#' @param \dots some methods for this generic require additional arguments.
#' None are used in this method.
#'
#' @export
summary.LMMsolve <- function(object,
                             which = "dimensions",
                             ...) {
  ## Checks.
  which <- match.arg(which)
  ## start and end of each variance component.
  Nres <- object$Nres
  e <- cumsum(object$varPar)
  s <- e - object$varPar + 1
  e <- e + Nres
  s <- s + Nres
  ## Get number of variance components.
  nVarComp <- length(object$term.labels.r)
  ## Get names of effective dimensions.
  namesED <- names(object$ED)
  ## Print total effective dimensions per variance component.
  for (i in 1:nVarComp) {
    ndx <- s[i]:e[i]
    cat(object$term.labels.r[i], "with total effective dimension ",
        round(sum(object$ED[ndx]), 2), "\n")
    if (length(ndx) > 1) {
      for (k in ndx) {
        cat("   ", namesED[k], "\t", round(object$ED[k], 2), "\n")
      }
    }
  }
}




