#' Fitted LMMsolve Object
#'
#' An object of class \code{LMMsolve} returned by the LMMsolve function,
#' representing a fitted linear mixed model. Objects of this class have
#' methods for the generic functions coef, fitted, residuals, loglik and
#' deviance.
#'
#' @return
#' An object of class \code{LMMsolve} contains the following components:
#' \item{logL}{The restricted log-likelihood at convergence}
#' \item{sigma2e}{The residual error}
#' \item{tau2e}{The estimated variance components}
#' \item{EDdf}{The effective dimensions}
#' \item{varPar}{The number of variance parameters for each variance component}
#' \item{VarDf}{The table with variance components}
#' \item{theta}{The precision parameters}
#' \item{coefMME}{A vector with all the estimated effects from mixed model equations}
#' \item{ndxCoefficients}{The indices of the coefficients with the names}
#' \item{yhat}{The fitted values}
#' \item{residuals}{The residuals}
#' \item{nIter}{The number of iterations for the mixed model to converge}
#' \item{y}{Response variable}
#' \item{X}{The design matrix for the fixed part of the mixed model}
#' \item{Z}{The design matrix for the random part of the mixed model}
#' \item{lGinv}{List with precision matrices for the random terms}
#' \item{lRinv}{List with precision matrices for the residual}
#' \item{C}{The mixed model coefficient matrix after last iteration}
#' \item{cholC}{The cholesky decomposition of coefficient matrix C}
#' \item{constantREML}{The REML constant}
#' \item{dim}{The dimensions for each of the fixed and random terms in the
#' mixed model}
#' \item{term.labels.f}{The names of the fixed terms in the mixed model}
#' \item{term.labels.r}{The names of the random terms in the mixed model}
#' \item{splRes}{An object with definition of spline argument}
#' \item{family}{An object of class family specifying the distribution and link function}
#' \item{trace}{A data.frame with the convergence sequence for the log likelihood and effective dimensions}.
#'
#' @usage NULL
#'
#' @rdname LMMsolveObject
#' @export
LMMsolveObject <- function(logL,
                           sigma2e,
                           tau2e,
                           EDdf,
                           varPar,
                           VarDf,
                           theta,
                           coefMME,
                           ndxCoefficients,
                           yhat,
                           residuals,
                           nIter,
                           y,
                           X,
                           Z,
                           lGinv,
                           lRinv,
                           C,
                           cholC,
                           constantREML,
                           dim,
                           Nres,
                           term.labels.f,
                           term.labels.r,
                           splRes,
                           family,
                           trace) {
  structure(list(logL = logL,
                 sigma2e = sigma2e,
                 tau2e = tau2e,
                 EDdf = EDdf,
                 varPar = varPar,
                 VarDf = VarDf,
                 theta = theta,
                 coefMME = coefMME,
                 ndxCoefficients = ndxCoefficients,
                 yhat = yhat,
                 residuals = residuals,
                 nIter = nIter,
                 y = y,
                 X = X,
                 Z = Z,
                 lGinv = lGinv,
                 lRinv = lRinv,
                 cholC = cholC,
                 C = C,
                 constantREML = constantREML,
                 dim = dim,
                 Nres = Nres,
                 term.labels.f = term.labels.f,
                 term.labels.r = term.labels.r,
                 splRes = splRes,
                 family = family,
                 trace = trace),
            class = c("LMMsolve", "list"))
}

#' Summarize Linear Mixed Model fits
#'
#' Summary method for class "LMMsolve". Creates either a table of effective
#' dimensions (which = "dimensions") or a table of variances (which =
#' "variances").
#'
#' @param object An object of class LMMsolve
#' @param which A character string indicating which summary table should be
#' created.
#' @param \dots Some methods for this generic require additional arguments.
#' None are used in this method.
#'
#' @return A data.frame with either effective dimensions or variances depending
#' on which.
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Obtain table of effective dimensions.
#' summ1 <- summary(LMM1)
#' print(summ1)
#'
#' ## Obtain table of variances.
#' summ2 <- summary(LMM1,
#'                 which = "variances")
#' print(summ2)
#'
#' @export
summary.LMMsolve <- function(object,
                             which = c("dimensions", "variances"),
                             ...) {
  ## Checks.
  which <- match.arg(which)
  if (which == "dimensions") {
    tbl <- object$EDdf
  } else if (which == "variances") {
    tbl <- object$VarDf
  }
  res <- structure(tbl,
                   class = c("summary.LMMsolve", "data.frame"),
                   which = which)
  return(res)
}

#' @param x An object of class summary.LMMsolve, the result of a call to
#' summary.LMM
#'
#' @describeIn summary.LMMsolve print summary
#'
#' @export
print.summary.LMMsolve <- function(x,
                                   ...) {
  which <- attr(x, which = "which")
  ## Compute sum of effective dimensions before rounding.
  EDTot <- sum(x[["Effective"]])
  ## Print max 2 decimals.
  x[2:ncol(x)] <- round(x[2:ncol(x)], 2)
  if (which == "dimensions") {
    cat("Table with effective dimensions and penalties: \n\n")
    print.data.frame(x, row.names = FALSE)
    cat("\n", "Total Effective Dimension:", EDTot, "\n")
  } else if (which == "variances") {
    cat("Table with variances: \n\n")
    print.data.frame(x, row.names = FALSE)
    cat("\n")
  }
}

#' Coefficients from the mixed model equations of an LMMsolve object.
#'
#' Obtain the coefficients from the mixed model equations of an LMMsolve object.
#'
#' @param object an object of class LMMsolve
#' @param se calculate standard errors, default FALSE.
#' @param \dots some methods for this generic require additional arguments.
#' None are used in this method.
#'
#' @return A list of vectors, containing the estimated effects for each fixed
#' effect and the predictions for each random effect in the defined linear
#' mixed model.
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Obtain coefficients.
#' coefs1 <- coef(LMM1)
#'
#' ## Obtain coefficients with standard errors.
#' coefs2 <- coef(LMM1, se = TRUE)
#'
#' @importFrom stats coef
#'
#' @export
coef.LMMsolve <- function(object,
                          se = FALSE,
                          ...) {
  u <- object$coefMME
  cf <- object$ndxCoefficients
  ## if not standard errors.
  if (!se) {
    coef <- lapply(X = cf, FUN = function(x) {
      ndx <- x != 0
      x[ndx] <- u[x[ndx]]
      return(x)
    })
  } else {
    ## if standard errors.
    n <- length(u)
    se <- calcStandardErrors(C = object$C,
                             D = spam::diag.spam(x = 1, nrow = n))
    coef <- lapply(X = cf, FUN = function(x) {
      x_se <- rep(NA, times = length(x))
      ndx <- x != 0
      xNdx <- x[ndx]
      x[ndx] <- u[xNdx]
      x_se[ndx] <- se[xNdx]
      df <- data.frame(coef = names(x), value = x, se = x_se,
                       zRatio = x / x_se, row.names = NULL)
      return(df)
    })
  }
  return(coef)
}

#' Fitted values of an LMMsolve object.
#'
#' Obtain the fitted values from a mixed model fitted using LMMSolve.
#'
#' @inheritParams coef.LMMsolve
#'
#' @return A vector of fitted values.
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Obtain fitted values.
#' fitted1 <- fitted(LMM1)
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
#' @return A vector of residuals.
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Obtain fitted values.
#' residuals1 <- residuals(LMM1)
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
#' @inheritParams coef.LMMsolve
#' @param includeConstant Should the constant in the restricted log-likelihood
#' be included. Default is \code{TRUE}, as for example in \code{lme4} and SAS.
#' In \code{asreml} the constant is omitted.
#'
#' @return The restricted maximum log-likelihood of the fitted model.
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Obtain log-likelihood.
#' logLik(LMM1)
#'
#' ## Obtain log-likelihood without constant.
#' logLik(LMM1, includeConstant = FALSE)
#'
#' @importFrom stats logLik
#'
#' @export
logLik.LMMsolve <- function(object,
                            includeConstant = TRUE,
                            ...) {
  logL <- object$logL
  if (includeConstant) {
    logL <- logL + object$constantREML
  }
  return(logL)
}

#' Deviance of an LMMsolve object
#'
#' Obtain the deviance of a model fitted using LMMsolve.
#'
#' @inheritParams logLik.LMMsolve
#'
#' @return The deviance of the fitted model.
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Obtain deviance.
#' logLik(LMM1)
#'
#' ## Obtain deviance. without constant.
#' logLik(LMM1, includeConstant = FALSE)
#'
#' @export
deviance.LMMsolve <- function(object,
                              includeConstant = TRUE,
                              ...) {
  logL <- logLik(object, includeConstant = includeConstant)
  dev <- -2 * logL
  return(dev)
}

#' Display the sparseness of the mixed model coefficient matrix
#'
#' @param object an object of class LMMsolve.
#' @param cholesky Should the cholesky decomposition of the coefficient matrix
#' be plotted?
#'
#' @return A plot of the sparseness of the mixed model coefficient matrix.
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Obtain deviance.
#' displayMME(LMM1)
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
#' @param object an object of class LMMsolve.
#'
#' @return A summary of the mixed model coefficient matrix and its choleski
#' decomposition.
#'
#' @examples
#' ## Fit model on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Obtain deviance.
#' diagnosticsMME(LMM1)
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



