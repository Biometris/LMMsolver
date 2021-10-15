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
  e <- cumsum(object$varPar) + Nres
  s <- e - object$varPar + 1
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

#' Coefficients from the mixed model equations of an LMMsolve object.
#'
#' Obtain the coefficients from the mixed model equations of an LMMsolve object.
#'
#' @param object an object of class LMMsolve
#' @param \dots some methods for this generic require additional arguments.
#' None are used in this method.
#'
#' @return A list of vectors, containing the estimated effects for each fixed
#' effect and the predictions for each random effect in the defined linear
#' mixed model.
#'
#' @export
coef.LMMsolve <- function(object,
                          ...) {
  result <- list()
  dim <- object$dim
  e <- cumsum(dim)
  s <- e - dim + 1

  for (i in 1:length(dim)) {
    result[[i]] <- object$a[s[i]:e[i]]
  }
  names(result) <- object$term.labels
  return(result)
}

#' Fitted values of an LMMsolve object.
#'
#' Obtain the fitted values from a mixed model fitted using LMMSolve.
#'
#' @inheritParams coef.LMMsolve
#'
#' @return A vector of fitted values.
#'
#' @export
fitted.LMMsolve <- function(object,
                            ...) {
  return(as.vector(object$yhat))
}

#' Residuals of an LMMsolve object.
#'
#' Obtain the residuals from a mixed model fitted using LMMSolve.
#'
#' @inheritParams coef.LMMsolve
#'
#' @return A vector of fitted values.
#'
#' @export
residuals.LMMsolve <- function(object,
                            ...) {
  return(as.vector(object$residuals))
}


#' Log-likelihood of an LMMsolve object
#'
#' Obtain the Restricted Maximum Log-Likelihood of a model fitted using
#' LMMsolve.
#'
#' @return The restricted maximum log-likelihood of the fitted model.
#'
#' @inheritParams coef.LMMsolve
#'
#' @export
logLik.LMMsolve <- function(object,
                            ...) {
  return(object$logL)
}

#' Deviance of an LMMsolve object
#'
#' Obtain the deviance of a model fitted using LMMsolve.
#'
#' @inheritParams coef.LMMsolve
#'
#' @return The deviance of the fitted model.
#'
#' @export
deviance.LMMsolve <- function(object,
                              ...) {
  return(object$dev)
}

#' Display the sparseness of the mixed model coefficient matrix
#'
#' @param object an object of class LMMsolve
#' @param cholesky logical. If \code{cholesky = TRUE} it will plot the cholesky.
#'
#' @export
displayMME <- function(object,
                       cholesky = FALSE) {
  if (!inherits(object, "LMMsolve")) {
    stop("object should be an object of class LMMsolve.\n")
  }
  if (!cholesky) {
    spam::display(object$C)
  } else {
    cholC <- chol(object$C)
    L <- t(spam::as.spam(cholC))
    spam::display(L)
  }
}

#' Give diagnostics for mixed model coefficient matrix C and the cholesky
#' decomposition
#'
#' @param object an object of class LMMsolve
#'
#' @export
diagnosticsMME <- function(object) {
  if (!inherits(object, "LMMsolve")) {
    stop("object should be an object of class LMMsolve.\n")
  }
  cat("Summary of matrix C \n")
  print(spam::summary.spam(object$C))
  cat("\n Summary of cholesky decomposition of C \n")
  print(spam::summary.spam(chol(object$C)))
}



