#' Solve Linear Mixed Models
#'
#' Solve Linear Mixed Models
#'
#' @param fixed ...
#' @param random ...
#' @param group ...
#' @param lGinverse ...
#' @param data ...
#' @param residual ...
#' @param tolerance ...
#' @param trace ...
#' @param display ...
#' @param maxit ...
#'
#' @importFrom stats model.frame terms model.matrix contrasts
#'
#' @export
LMMsolve <- function(fixed,
                     random = NULL,
                     group = NULL,
                     lGinverse = NULL,
                     data,
                     residual = NULL,
                     tolerance = 1.0e-6,
                     trace = FALSE,
                     display = FALSE,
                     maxit = 250) {
  ## make random part:
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
    Z2 <- data[, ndx]
  } else {
    dim2.r <- NULL
    term2.labels.r <- NULL
    Z2 <- NULL
  }

  if (!(is.null(random) & is.null(group))) {
    if (is.null(random)) {
      Z <- Z2
    } else if (is.null(group)) {
      Z <- Z1
    } else {
      Z <- cbind(Z1, Z2)
    }
    Z <- as.matrix(Z)
    dim.r <- c(dim1.r, dim2.r)
    term.labels.r <- c(term1.labels.r, term2.labels.r)

    e <- cumsum(dim.r)
    s <- e - dim.r + 1

    lGinv <- list()
    for(i in 1:length(dim.r)) {
      if (term.labels.r[i] %in% names(lGinverse)) {
        #print(term.labels.r[i])
        #lGinv[[i]] <-
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

  ## make fixed part:
  mf <- model.frame(fixed, data)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
  X = model.matrix(mt, data = mf,
                   contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                          FUN = contrasts, contrasts = TRUE))
  dim.f <- table(attr(X, "assign"))
  term.labels.f <- attr(mt, "term.labels")

  # add intercept....
  if (attr(mt, "intercept") == 1) {
    term.labels.f <- c("(Intercept)", term.labels.f)
  }

  if (!is.null(residual)) {
    lRinv <- makeRlist(df = data, column = residual)
  } else {
    lRinv <- list()
    n <- nrow(data)
    lRinv[[1]] <- spam::diag.spam(1, n)
    names(lRinv) = "residual"
  }
  y <- mf[, 1]
  obj <- sparseMixedModels(y, X, Z, lGinv, lRinv, tolerance, trace,
                           display = display, maxit = maxit)
  dim <- as.numeric(c(dim.f, dim.r))
  term.labels <- c(term.labels.f, term.labels.r)
  obj$dim <- dim
  obj$term.labels <- term.labels
  class(obj) <- "LMMsolve"
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
