#' Solve Linear Mixed Models
#'
#' Solve Linear Mixed Models
#'
#' @param fixed A formula for the fixed part of the model. Should be of the
#' form "response ~ pred"
#' @param random A formula for the random part of the model. Should be of the
#' form "~ pred".
#' @param spatial A formula for the spatial part of the model. Should be of the
#' form "~ sap2D()" or "~sap3d()".
#' @param group A named list where each component is a numeric vector
#' specifying contiguous fields in data that are to be considered as a
#' single term.
#' @param lGinverse A named list where each component is a matrix corresponding
#' to the group of the same name.
#' @param data A data.frame containing the modeling data.
#' @param residual A formula for the residual part of the model. Should be of
#' the form "~ pred".
#' @param tolerance A numerical value. The convergence tolerance for the
#' algorithm.
#' @param trace Should the progress of the algorithm be printed?
#' @param display Should the ... matrix created in the algorithm be plotted?
#' @param maxit A numerical value. The maximum number of iterations for the
#' algorithm.
#'
#' @return An object of class LMMsolve, a list with the following items:
#' \item{logL}{The loglikelihood}
#' \item{sigma2e}{The residual error}
#' \item{tau2e}{}
#' \item{ED}{The effective dimensions}
#' \item{EDmax}{The maximal effective dimensions}
#' \item{EDnames}{The names of the effective dimensions}
#' \item{a}{}
#' \item{yhat}{The fitted values}
#' \item{dim}{}
#' \item{term.labels}{}
#'
#' @importFrom stats model.frame terms model.matrix contrasts as.formula terms.formula
#'
#' @export
LMMsolve <- function(fixed,
                     random = NULL,
                     spatial = NULL,
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
  if (!is.null(spatial) && length(terms(spatial)) != 2) {
    stop("spatial model formula must be of form \"~ sap2D()\" or \"~sap3d()\".\n")
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
  mf <- model.frame(fixed, data)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
  X = model.matrix(mt, data = mf,
                   contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                          FUN = contrasts, contrasts = TRUE))
  dim.f <- table(attr(X, "assign"))
  term.labels.f <- attr(mt, "term.labels")

  ## Make spatial part.
  if (!is.null(spatial)) {

    if (inherits(spatial, "character")) {
      spatial <- as.formula(spatial)
    }
    tf <- terms(spatial, specials = c("sap2D", "sap3D"))
    terms <- attr(tf, "term.labels")
    nt <- length(terms)

    spatRes <- eval(parse(text = terms), envir = data, enclos = parent.frame())

    X <- cbind(X, spatRes$X)
    Z <- cbind(Z, spatRes$Z)
    lGinv <- c(lGinv, spatRes$lGinv)
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
  class(obj) <- c("LMMsolve", "list")
  return(obj)
}

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
