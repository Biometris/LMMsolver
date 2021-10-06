#' Solve Linear Mixed Models
#'
#' Solve Linear Mixed Models.
#'
#' @param fixed A formula for the fixed part of the model. Should be of the
#' form "response ~ pred"
#' @param random A formula for the random part of the model. Should be of the
#' form "~ pred".
#' @param spline A formula for the spline part of the model. Should be of the
#' form "~ spl1D()", ~ spl2D()" or "~spl3D()".
#' @param group A named list where each component is a numeric vector
#' specifying contiguous fields in data that are to be considered as a
#' single term.
#' @param data A data.frame containing the modeling data.
#' @param residual A formula for the residual part of the model. Should be of
#' the form "~ pred".
#' @param tolerance A numerical value. The convergence tolerance for the
#' modified Henderson algorithm to estimate the variance components.
#' @param trace Should the progress of the algorithm be printed? Default
#' \code{trace = FALSE}.
#' @param maxit A numerical value. The maximum number of iterations for the
#' algorithm. Default \code{maxit = 250}.
#' @param omitConstant Omit the constant in the restricted log-likelihood. Default is
#' \code{TRUE}, as for example in \code{asreml}. In \code{nlme} and SAS the constant
#' is included.
#'
#' @return An object of class \code{LMMsolve} representing the fitted model.
#' See \code{\link{LMMsolveObject}} for a full description of the components in
#' this object.
#'
#' @seealso \code{\link{LMMsolveObject}}
#'
#' @importFrom stats model.frame terms model.matrix contrasts as.formula
#' terms.formula
#'
#' @export
LMMsolve <- function(fixed,
                     random = NULL,
                     spline = NULL,
                     group = NULL,
                     data,
                     residual = NULL,
                     tolerance = 1.0e-6,
                     trace = FALSE,
                     maxit = 250,
                     omitConstant = TRUE) {
  ## Input checks.
  if (!inherits(data, "data.frame")) {
    stop("data should be a data.frame.\n")
  }
  if (!inherits(fixed, "formula") || length(terms(fixed)) != 3) {
    stop("fixed should be a formula of the form \"resp ~ pred\".\n")
  }
  if (!is.null(random) &&
      (!inherits(random, "formula") || length(terms(random)) != 2)) {
    stop("random should be a formula of the form \" ~ pred\".\n")
  }
  if (!is.null(spline) &&
      (!inherits(spline, "formula") || length(terms(spline)) != 2 ||
       ## Spline formula should consist of splxD() and nothing else.
       sum(!sapply(attr(terms(spline, specials = c("spl1D", "spl2D", "spl3D")),
                        "specials"), is.null)) != 1)) {
    stop("spline should be a formula of form \"~ spl1D()\", \"~ spl2D()\" ",
         "or \"~spl3D()\".\n")
  }
  if (!is.null(residual) &&
      (!inherits(residual, "formula") || length(terms(residual)) != 2)) {
    stop("residual should be a formula of the form \" ~ pred\".\n")
  }
  if (!is.numeric(tolerance) || length(tolerance) > 1 || tolerance < 0) {
    stop("tolerance should be a positive numerical value.")
  }
  if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 0) {
    stop("maxit should be a positive numerical value.")
  }
  ## Check that all variables used in formulas are in data.
  checkFormVars(fixed, data)
  checkFormVars(random, data)
  checkFormVars(residual, data)

  ## Remove NA for response variable from data.
  respVar <- all.vars(fixed)[attr(terms(fixed), "response")]
  respVarNA <- is.na(data[[respVar]])
  if (sum(respVarNA) > 0) {
    warning(sum(respVarNA), " observations removed with missing value for ",
            respVar, ".\n", call. = FALSE)
    data <- data[!respVarNA, ]
  }

  ## Make random part.
  if (!is.null(random)) {
    mf <- model.frame(random, data, drop.unused.levels = TRUE, na.action = NULL)
    mt <- terms(mf)
    f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
    Z1 <- model.matrix(mt, data = mf,
                       contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                              FUN = contrasts,
                                              contrasts = FALSE))
    dim1.r <- table(attr(Z1, "assign"))[-1]
    term1.labels.r <- attr(mt, "term.labels")
    ## Number of variance parameters (see Gilmour 1995) for each variance component
    varPar1 <- rep(1, length(dim1.r))
    Z1 <- Z1[, -1]
  } else {
    dim1.r <- NULL
    term1.labels.r <- NULL
    Z1 <- NULL
    varPar1 <- NULL
  }

  if (!is.null(group)) {
    ndx <- unlist(group)
    dim2.r <- sapply(X = group, FUN = length)
    term2.labels.r <- names(group)
    varPar2 <- rep(1, length(dim2.r))
    Z2 <- as.matrix(data[, ndx])
  } else {
    dim2.r <- NULL
    term2.labels.r <- NULL
    Z2 <- NULL
    varPar2 <- NULL
  }

  if (!(is.null(random) & is.null(group))) {
    Z <- cbind(Z1, Z2)
    dim.r <- c(dim1.r, dim2.r)
    term.labels.r <- c(term1.labels.r, term2.labels.r)
    varPar <- c(varPar1, varPar2)
    e <- cumsum(dim.r)
    s <- e - dim.r + 1

    lGinv <- list()
    for (i in 1:length(dim.r)) {
      tmp <- rep(0, sum(dim.r))
      tmp[s[i]:e[i]] <- 1
      lGinv[[i]] <- spam::diag.spam(tmp)
    }
    names(lGinv) <- term.labels.r
  } else {
    Z <- NULL
    lGinv <- NULL
    dim.r <- NULL
    term.labels.r <- NULL
    varPar <- NULL
  }

  ## Make fixed part.
  mf <- model.frame(fixed, data, drop.unused.levels = TRUE)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
  X = model.matrix(mt, data = mf,
                   contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                          FUN = contrasts, contrasts = TRUE))
  dim.f <- table(attr(X, "assign"))
  term.labels.f <- attr(mt, "term.labels")

  ## Add spline part.
  splRes <- NULL
  if (!is.null(spline)) {
    tf <- terms(spline, specials = c("spl1D", "spl2D", "spl3D"))
    terms <- attr(tf, "term.labels")
    splRes <- eval(parse(text = terms), envir = data, enclos = parent.frame())
    ## Add to design matrix fixed effect X.
    X <- cbind(X, splRes$X)
    ## Add to design matrix random effect Z.
    Z <- cbind(Z, splRes$Z)
    ## Expand matrices Ginv to the updated Z.
    lGinv <- expandGinv(lGinv, splRes$lGinv)
    ## A splxD model has x parameters.
    varPar <- c(varPar, length(splRes$lGinv))
    ## Add dims.
    dim.f <- c(dim.f, splRes$dim.f)
    dim.r <- c(dim.r, splRes$dim.r)
    ## Add labels.
    term.labels.f <- c(term.labels.f, splRes$term.labels.f)
    term.labels.r <- c(term.labels.r, splRes$term.labels.r)
  }

  ## Add intercept.
  if (attr(mt, "intercept") == 1) {
    term.labels.f <- c("(Intercept)", term.labels.f)
  }

  ## Make residual part.
  if (!is.null(residual)) {
    residVar <- all.vars(residual)
    lRinv <- makeRlist(df = data, column = residVar)
  } else {
    lRinv <- list(residual = spam::diag.spam(1, nrow(data)))
    attr(lRinv, "cnt") <- nrow(data)
  }
  y <- mf[, 1]
  obj <- sparseMixedModels(y = y, X = X, Z = Z, lGinv = lGinv, lRinv = lRinv,
                           tolerance = tolerance, trace = trace, maxit = maxit)
  NomEffDimRes <- attr(lRinv,"cnt") - 1
  NomEffDimRan <- calcNomEffDim(X, Z, dim.r)
  NomEffDim <- c(NomEffDimRes, NomEffDimRan)

  # not working: names(NomEffDim) <- names(obj$ED),
  # as NomEffDim is per variance component:
  obj$EDnominal <- NomEffDim

  if (!omitConstant)
  {
    Nobs <- length(y)
    p <- sum(dim.f)
    Constant = -0.5*log(2*pi)*(Nobs-p)
    obj$logL <- obj$logL + Constant
  }
  dim <- as.numeric(c(dim.f, dim.r))
  term.labels <- c(term.labels.f, term.labels.r)
  obj$varPar <- varPar
  obj$dim <- dim
  obj$Nres <- length(lRinv)
  obj$term.labels.f <- term.labels.f
  obj$term.labels.r <- term.labels.r
  obj$term.labels <- term.labels
  obj$splRes <- splRes
  obj$dev <- -2.0*obj$logL
  return(LMMsolveObject(obj))
}

#' Obtain the coefficients from the mixed model equations
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

  for(i in 1:length(dim)) {
    result[[i]] <- object$a[s[i]:e[i]]
  }
  names(result) <- object$term.labels
  return(result)
}


#' Display the sparseness of the mixed model coefficient matrix.
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
#' decomposition.
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


