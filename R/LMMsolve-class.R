#' Fitted LMMsolve Object
#'
#' An object of class \code{LMMsolve} returned by the LMMsolve function,
#' representing a fitted linear mixed model. Objects of this class have
#' methods for the generic functions coef, fitted, residuals, loglik and
#' deviance.
#'
#' @returns
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
#' \item{respVar}{The name(s) of the response variable(s).}
#' \item{splRes}{An object with definition of spline argument}
#' \item{deviance}{The relative deviance}
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
                           respVar,
                           splRes,
                           family,
                           deviance,
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
                 respVar = respVar,
                 splRes = splRes,
                 family = family,
                 deviance = deviance,
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
#' @returns A data.frame with either effective dimensions or variances depending
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
#' @returns A list of vectors, containing the estimated effects for each fixed
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

  ## remove termType attribute, only needed for predict function
  cf <- lapply(cf, function(x) { attr(x,which="termType") <- NULL
                                 return(x) })

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
#' @returns A vector of fitted values.
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
#' @returns A vector of residuals.
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
#' @returns The restricted maximum log-likelihood of the fitted model.
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
#' @param relative Deviance relative conditional or absolute unconditional
#' (-2*logLik(object))? Default \code{relative = TRUE}.
#'
#' @returns The deviance of the fitted model.
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
#' deviance(LMM1)
#'
#' @export
deviance.LMMsolve <- function(object,
                              relative = TRUE,
                              includeConstant = TRUE,
                              ...) {
  if (relative) {
    return(object$deviance)
  }
  # else
  logL <- logLik(object, includeConstant = includeConstant)
  dev <- -2 * logL
  return(dev)
}

#' Predict function
#'
#' @param object an object of class LMMsolve.
#' @param ... Unused.
#' @param newdata A data.frame containing new points for which the smooth
#' trend should be computed. Column names should include the names used when
#' fitting the spline model.
#' @param se.fit calculate standard errors, default \code{FALSE}.
#'
#' @returns A data.frame with predictions for the smooth trend on the specified
#' grid. The standard errors are saved if `se.fit=TRUE`.
#'
#' @examples
#' ## simulate some data
#' f <- function(x) { 0.3 + 0.4*x + 0.2*sin(20*x) }
#' set.seed(12)
#' n <- 150
#' x <- seq(0, 1, length = n)
#' sigma2e <- 0.04
#' y <- f(x) + rnorm(n, sd = sqrt(sigma2e))
#' dat <- data.frame(x, y)
#'
#' ## fit the model
#' obj <- LMMsolve(fixed = y ~ 1,
#'          spline = ~spl1D(x, nseg = 50), data = dat)
#'
#' ## make predictions on a grid
#' newdat <- data.frame(x = seq(0, 1, length = 300))
#' pred <- predict(obj, newdata = newdat, se.fit = TRUE)
#' head(pred)
#'
#' @export
predict.LMMsolve <- function(object,
                             ...,
                             newdata,
                             se.fit = FALSE) {
  if (!inherits(object, "LMMsolve")) {
    stop("object should be an object of class LMMsolve.\n")
  }
  #if (is.null(object$splRes)) {
  #  stop("The model was fitted without a spline component.\n")
  #}
  if (!inherits(newdata, "data.frame")) {
    stop("newdata should be a data.frame.\n")
  }
  family <- object$family
  if (family$family == "multinomial") {
    if (se.fit) {
      stop("se.fit=TRUE not implemented yet for multinomial.\n")
    }
    ndxCf <- object$ndxCoefficients
    IsFactor <- sapply(ndxCf, FUN=function(x) {return(attr(x,"termType")=="factor")})
    if (any(IsFactor)) {
      stop("use of factors not implemented yet for multinomial.\n")
    }
  }

  varNames <- unlist(sapply(object$splRes,function(z){names(z$x)}))
  colNames <- colnames(newdata)
  Missing <- !(varNames %in% colNames)
  if (sum(Missing) > 0) {
    missingVar <- paste(varNames[Missing], collapse=",")
    str <- paste0("variables (", missingVar, ") in data.frame newdata missing.\n")
    stop(str)
  }

  # some items not implemented yet
  #s1 <-sum(sapply(obj2$ndxCoefficients,
  #                 function(x) {attr(x,which="termType") == "factor"}))
  #if (s1 > 0) stop("predict function for factors not implemented yet")
  s2 <-sum(sapply(object$ndxCoefficients,
                   function(x) {attr(x,which="termType") == "grp"}))
  if (s2 > 0) stop("predict function for grp() not implemented yet")

  splFixLab <- sapply(object$splRes, function(x) { x$term.labels.f })
  splRanLab <- sapply(object$splRes, function(x) { x$term.labels.r })

  nGam <- length(object$splRes)
  xGrid <- list()
  BxTot <- list()
  XTot <- list()
  for (s in seq_len(nGam)) {
    spl <- object$splRes[[s]]
    ## check whether the values are in range Bsplines
    chkValBsplines(spl, newdata)
    x <- spl$x
    xGrid[[s]] <- lapply(X = seq_along(x), FUN = function(i) {
      newdata[[names(x)[i]]]})
    Bx <- mapply(FUN = Bsplines, spl$knots, xGrid[[s]], deriv=0)
    BxTot[[s]] <- Reduce(RowKronecker, Bx)
    G <- lapply(X=spl$knots, FUN = function(x) {
      constructG(knots = x, scaleX = spl$scaleX, pord = spl$pord)})
    ## Compute G over all dimensions
    GTot <- Reduce('%x%', G)
    ## no scaling for first column of GTot
    GTot[,1] <- 1
    X <- BxTot[[s]] %*% GTot
    ## Remove intercept (needed when fitting model to avoid singularities).
    XTot[[s]] <- removeIntercept(X)
  }
  dim <- object$dim
  lU <- list()
  nRow <- nrow(newdata)
  if (family$family == "multinomial") {
    nCat <- length(object$respVar) - 1
  } else {
    nCat <- 1
  }

  for (i in seq_along(dim)) {
    lU[[i]] = spam::spam(x = 0, nrow = nRow, ncol = dim[i]/nCat)
  }
  # intercept:
  lU[[1]] = spam::spam(x = 1, nrow = nRow, ncol = 1)

  labels <- c(object$term.labels.f, object$term.labels.r)

  for (s in seq_len(nGam)) {
    spl <- object$splRes[[s]]
    ndx.f <- which(spl$term.labels.f == labels)
    ndx.r <- which(spl$term.labels.r == labels)
    lU[[ndx.f]] <- XTot[[s]]
    lU[[ndx.r]] <- BxTot[[s]]
  }
  U <- Reduce(spam::cbind.spam, lU)
  if (family$family == "multinomial") {
    U <- U %x% spam::diag.spam(nCat)
  }

  tmp <- object$term.labels.f[-1]
  fixTerms <- setdiff(tmp, splFixLab)

  nFixTerms <- length(fixTerms)
  if (nFixTerms > 0) {
    colNames <- colnames(newdata)
    Missing <- !(fixTerms %in% colNames)
    if (sum(Missing) > 0) {
      missingVar <- paste(fixTerms[Missing], collapse=",")
      str <- paste0("variables (", missingVar, ") in data.frame newdata missing.\n")
      stop(str)
    }
    for (i in seq_len(nFixTerms)) {
      U <- U + makeDesignTerm(object, newdata, fixTerms[i])
    }
  }

  outDat <- newdata

  ranTerms <- setdiff(object$term.labels.r, splRanLab)
  nRanTerms <- length(ranTerms)
  for (i in seq_len(nRanTerms)) {
    term <- ranTerms[[i]]
    outDat[[term]] <- rep("Excluded",nRow)
  }
  if (family$family == "multinomial") {
    eta0 <- as.vector(U %*% object$coefMME)
    etaM <- matrix(data=eta0,nrow = nRow, ncol= nCat, byrow=TRUE)
    pi_est <- t(apply(etaM, MARGIN=1, FUN = family$linkinv)) # linkinv
    pi_est <- cbind(pi_est, 1.0 - rowSums(pi_est))
    colnames(pi_est) <- object$respVar
    tmp <- NULL
    for (i in seq_along(object$respVar)) {
      tmp2 <- data.frame(outDat, category = object$respVar[i])
      tmp <- rbind(tmp, tmp2)
    }
    outDat <- tmp
    outDat[["category"]] <- as.factor(rep(object$respVar, each = nRow))
    familyPred <- as.vector(pi_est)
  } else {
    eta <- as.vector(U %*% object$coefMME)
    family <- object$family
    familyPred <- family$linkinv(eta)
  }

  outDat[["ypred"]] <- familyPred
  if (se.fit) {
    outDat[["se"]] <- calcStandardErrors(C = object$C, D = U)*abs(family$mu.eta(eta))
  }

  return(outDat)
}




