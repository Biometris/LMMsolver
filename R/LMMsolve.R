#' Solve Linear Mixed Models
#'
#' Solve Linear Mixed Models
#'
#' @param fixed A formula for the fixed part of the model. Should be of the
#' form "response ~ pred"
#' @param random A formula for the random part of the model. Should be of the
#' form "~ pred".
#' @param spline A formula for the spline part of the model. Should be of the
#' form "~ spl1D()", ~ sap2D()" or "~sap3d()".
#' @param group A named list where each component is a numeric vector
#' specifying contiguous fields in data that are to be considered as a
#' single term.
#' @param lGinverse A named list where each component is a matrix corresponding
#' to the group of the same name.
#' @param data A data.frame containing the modeling data.
#' @param residual A formula for the residual part of the model. Should be of
#' the form "~ pred".
#' @param tolerance A numerical value. The convergence tolerance for the
#' modified Henderson algorithm to estimate the variance components.
#' @param trace Should the progress of the algorithm be printed? Default \code{trace=FALSE}.
#' @param display Should the sparse matrix created in the algorithm be plotted?
#' @param maxit A numerical value. The maximum number of iterations for the
#' algorithm. Default \code{maxit=250}.
#'
#' @return An object of class LMMsolve, a list with the following items:
#' \item{logL}{The loglikelihood}
#' \item{sigma2e}{The residual error}
#' \item{tau2e}{estimated variance components}
#' \item{ED}{The effective dimensions}
#' \item{EDmax}{The maximal effective dimensions}
#' \item{EDnames}{The names of the effective dimensions}
#' \item{a}{the estimated effects from the mixed model equations}
#' \item{yhat}{The fitted values}
#' \item{dim}{dimensions for each fixed or random term in the mixed model}
#' \item{term.labels}{names of the fixed and random terms in the mixed model}
#' \item{splRes}{An object with definition of spline argument}
#' @importFrom stats model.frame terms model.matrix contrasts as.formula terms.formula
#'
#' @export
LMMsolve <- function(fixed,
                     random = NULL,
                     spline = NULL,
                     group = NULL,
                     lGinverse = NULL,
                     data,
                     residual = NULL,
                     tolerance = 1.0e-6,
                     trace = FALSE,
                     display = FALSE,
                     maxit = 250) {
  ## Input checks.
  if (!inherits(data, "data.frame")) {
    stop("data should be a data.frame.\n")
  }
  if (length(terms(fixed)) != 3) {
    stop("fixed model formula must be of the form \"resp ~ pred\".\n")
  }
  if (!is.null(random) && length(terms(random)) != 2) {
    stop("random model formula must be of form \" ~ pred\".\n")
  }
  if (!is.null(spline) && length(terms(spline)) != 2) {
    stop("spline model formula must be of form \"~ sap2D()\" or \"~sap3d()\".\n")
  }
  if (!is.null(residual) && length(terms(residual)) != 2) {
    stop("residual model formula must be of the form \" ~ pred\".\n")
  }
  if (!is.numeric(tolerance) || length(tolerance) > 1 || tolerance < 0) {
    stop("tolerance should be a positive numerical value.")
  }
  if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 0) {
    stop("maxit should be a positive numerical value.")
  }

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
    Z1 <- Z1[, -1]
  } else {
    dim1.r <- NULL
    term1.labels.r <- NULL
    Z1 <- NULL
  }

  if (!is.null(group)) {
    ndx <- unlist(group)
    dim2.r <- sapply(X = group, FUN = length)
    term2.labels.r <- names(group)
    Z2 <- as.matrix(data[, ndx])
  } else {
    dim2.r <- NULL
    term2.labels.r <- NULL
    Z2 <- NULL
  }

  if (!(is.null(random) & is.null(group))) {
    Z <- cbind(Z1, Z2)
    dim.r <- c(dim1.r, dim2.r)
    term.labels.r <- c(term1.labels.r, term2.labels.r)

    e <- cumsum(dim.r)
    s <- e - dim.r + 1

    lGinv <- list()
    for(i in 1:length(dim.r)) {
      if (term.labels.r[i] %in% names(lGinverse)) {
        tmp <- spam::diag.spam(0, sum(dim.r))
        tmp[s[i]:e[i], s[i]:e[i]] <- lGinverse[[term.labels.r[i]]]
        lGinv[[i]] <- tmp
      } else {
        tmp <- rep(0, sum(dim.r))
        tmp[s[i]:e[i]] <- 1
        lGinv[[i]] <- spam::diag.spam(tmp)
      }
    }
    names(lGinv) <- term.labels.r
  } else {
    Z <- NULL
    lGinv <- NULL
    dim.r <- NULL
    term.labels.r <- NULL
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

    if (inherits(spline, "character")) {
      spline <- as.formula(spline)
    }
    tf <- terms(spline, specials = c("spl1D", "sap2D", "sap3D"))
    terms <- attr(tf, "term.labels")
    nt <- length(terms)

    splRes <- eval(parse(text = terms), envir = data, enclos = parent.frame())

    ## Add to design matrix fixed effect X
    X <- cbind(X, splRes$X)
    ## Add to design matrix random effect Z
    Z <- cbind(Z, splRes$Z)
    ## Expand matrices Ginv to the updated Z
    lGinv <- ExpandGinv(lGinv, splRes$lGinv)

    ## Add dims
    dim.f <- c(dim.f, splRes$dim.f)
    dim.r <- c(dim.r, splRes$dim.r)
    ## Add labels
    term.labels.f <- c(term.labels.f, splRes$term.labels.f)
    term.labels.r <- c(term.labels.r, splRes$term.labels.r)

  }

  ## Add intercept.
  if (attr(mt, "intercept") == 1) {
    term.labels.f <- c("(Intercept)", term.labels.f)
  }

  ## Make residual part.
  if (!is.null(residual)) {
    lRinv <- makeRlist(df = data, column = residual)
  } else {
    lRinv <- list(residual = spam::diag.spam(1, nrow(data)))
  }
  y <- mf[, 1]
  obj <- sparseMixedModels(y = y, X = X, Z = Z, lGinv = lGinv, lRinv = lRinv,
                           tolerance = tolerance, trace = trace,
                           display = display, maxit = maxit)
  dim <- as.numeric(c(dim.f, dim.r))
  term.labels <- c(term.labels.f, term.labels.r)
  obj$dim <- dim
  obj$term.labels <- term.labels
  obj$splRes <- splRes
  class(obj) <- c("LMMsolve", "list")
  return(obj)
}

#' Obtain the coefficients from the mixed model equations
#'
#' @param object an object of class LMMsolve
#' @result a list of vectors, containing the estimated effects for each fixed effect
#' and the predictions for each random effect in the defined linear mixed model.
#'
#' @export
coef.LMMsolve <- function(object) {
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
