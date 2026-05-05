build_random_Z1_train <- function(random, data) {

  if (is.null(random)) {
    return(list(Z = NULL, spec = NULL))
  }

  mf <- model.frame(random, data, drop.unused.levels = TRUE)

  ## name check
  names_mf <- names(mf)
  if (!all(names_mf == make.names(names_mf))) {
    stop("Syntactically invalid name(s): ",
         paste(names_mf[names_mf != make.names(names_mf)], collapse=", "),
         call. = FALSE)
  }

  mt <- terms(mf)

  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]

  contrasts.arg <- lapply(
    mf[, f.terms, drop = FALSE],
    contrasts,
    contrasts = FALSE   # important for random
  )

  Z1 <- Matrix::sparse.model.matrix(
    mt,
    data = mf,
    contrasts.arg = contrasts.arg
  )

  dim.r <- table(attr(Z1, "assign"))[-1]
  term.labels.r <- attr(mt, "term.labels")

  ## remove intercept column
  if (ncol(Z1) > 1) {
    Z1 <- Z1[, -1, drop = FALSE]
    colnames_Z1 <- colnames(Z1)
    Z1 <- spam::as.spam.dgCMatrix(Z1)
  } else {
    Z1 <- NULL
    colnames_Z1 <- NULL
  }

  if (!is.null(Z1)) {
    spec <- list(
      terms        = mt,
      xlevels      = .getXlevels(mt, mf),
      contrasts    = contrasts.arg,
      colnames     = colnames_Z1,
      dim.r        = dim.r,
      term.labels  = term.labels.r)
  } else {
    spec <- list(
      terms       = NULL,
      xlevels     = NULL,
      constrasts  = NULL,
      colnames    = NULL,
      dim.r       = NULL,
      term.labels = NULL)
  }

  return(list(Z = Z1, spec = spec))
}

constructRandom <- function(random, group, condFactor, data) {
  res <- build_random_Z1_train(random, data)

  Z1 <- res$Z
  spec_Z1 <- res$spec
  dim1.r <- spec_Z1$dim.r
  term1.labels.r <- spec_Z1$term.labels
  if (!is.null(Z1)) {
    scFactor1 <- rep(1, length(dim1.r))
    ## Number of variance parameters (see Gilmour 1995) for each variance component
    varPar1 <- rep(1, length(dim1.r))
  } else {
    scFactor1 <- NULL
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
    for (i in seq_along(dim.r)) {
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
            varPar = varPar,
            nNonSplinesRandom = length(dim.r),
            ran.spec = spec_Z1)
  return(L)
}


build_random_Z1_pred <- function(spec, data) {

  if (is.null(spec)) return(NULL)

  mt <- spec$terms

  ## 1. Ensure required variables exist (inject NA if missing)
  data <- fill_missing_vars(mt, data, spec$xlevels)

  ## 2. Check for new (invalid) levels, ignoring NA
  check_new_levels(data, spec$xlevels)

  ## 3. Build model frame (keep NA)
  mf <- model.frame(
    mt,
    data,
    xlev = spec$xlevels,
    drop.unused.levels = FALSE,
    na.action = na.pass
  )

  ## 4. Construct design matrix
  Z1 <- Matrix::sparse.model.matrix(
    mt,
    data = mf,
    contrasts.arg = spec$contrasts
  )

  ## 5. Convert NA → 0 (missing = no contribution)
  if (anyNA(Z1)) {
    Z1[is.na(Z1)] <- 0
  }

  ## 6. Remove intercept column (as in training)
  if (ncol(Z1) > 1) {
    Z1 <- Z1[, -1, drop = FALSE]
  } else {
    return(NULL)
  }

  ## 7. make sure column order is correct
  Z1 <- Z1[, spec$colnames, drop = FALSE]

  ## 8. Convert to spam
  Z1 <- spam::as.spam.dgCMatrix(Z1)

  return(Z1)
}

