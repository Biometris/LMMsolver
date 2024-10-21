constructRandom <- function(random, group, condFactor, data) {
  ## Make random part.
  if (!is.null(random)) {
    mf <- model.frame(random, data, drop.unused.levels = TRUE, na.action = NULL)
    mt <- terms(mf)
    f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
    Z1 <- Matrix::sparse.model.matrix(mt, data = mf,
                                      contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                                             FUN = contrasts,
                                                             contrasts = FALSE))
    dim1.r <- table(attr(Z1, "assign"))[-1]
    term1.labels.r <- attr(mt, "term.labels")
    scFactor1 <- rep(1, length(dim1.r))
    ## Number of variance parameters (see Gilmour 1995) for each variance component
    varPar1 <- rep(1, length(dim1.r))
    if (ncol(Z1) > 1) {
      Z1 <- Z1[, -1, drop = FALSE]
      Z1 <- spam::as.spam.dgCMatrix(Z1)
    }
    else {
      Z1 <- NULL
    }
  } else {
    scFactor1 <- NULL
    dim1.r <- NULL
    term1.labels.r <- NULL
    Z1 <- NULL
    varPar1 <- NULL
  }
  if (!is.null(group)) {
    ndx <- unlist(group)
    dim2.r <- sapply(X = group, FUN = length)
    term2.labels.r <- names(group)
    scFactor2 <- rep(1, length(dim2.r))
    varPar2 <- rep(1, length(dim2.r))
    Z2 <- spam::as.spam(as.matrix(data[, ndx]))
  } else {
    dim2.r <- NULL
    term2.labels.r <- NULL
    scFactor2 <- NULL
    Z2 <- NULL
    varPar2 <- NULL
  }
  if (!is.null(condFactor)) {
    dim3.r <- condFactor$dim.r
    term3.labels.r <- condFactor$term.labels.r
    scFactor3 <- rep(1, length(dim3.r))
    varPar3 <- rep(1, length(dim3.r))
    Z3 <- condFactor$Z
  } else {
    dim3.r <- NULL
    term3.labels.r <- NULL
    scFactor3 <- NULL
    Z3 <- NULL
    varPar3 <- NULL
  }
  if (!(is.null(random) && is.null(group) && is.null(condFactor))) {
    Z <- spam::cbind.spam(Z1, Z2, Z3)
    dim.r <- c(dim1.r, dim2.r, dim3.r)
    term.labels.r <- c(term1.labels.r, term2.labels.r, term3.labels.r)
    scFactor <- c(scFactor1, scFactor2, scFactor3)
    varPar <- c(varPar1, varPar2, varPar3)
    e <- cumsum(dim.r)
    s <- e - dim.r + 1
    lGinv <- list()
    for (i in 1:length(dim.r)) {
      tmp <- rep(0, sum(dim.r))
      tmp[s[i]:e[i]] <- 1
      lGinv[[i]] <- spam::cleanup(spam::diag.spam(tmp))
    }
    names(lGinv) <- term.labels.r
  } else {
    Z <- NULL
    lGinv <- NULL
    dim.r <- NULL
    term.labels.r <- NULL
    scFactor <- NULL
    varPar <- NULL
  }
  L <- list(Z = Z, lGinv = lGinv, dim.r = dim.r,
            term.labels.r = term.labels.r, scFactor = scFactor,
            varPar = varPar)
  return(L)
}
