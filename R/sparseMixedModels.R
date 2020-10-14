
# mixed model coefficient matrix...

# calcalate the weighted sum of list of matrices.
linearSum  <- function(theta,
                       matrixList) {
  C <- Reduce('+', Map(f = function(x, y) x * y, theta, matrixList))
  return(C)
}

# calculate quadratic form xAx (efficient):
quadForm <- function(x,
                     A,
                     y = x) {
  sum(x*(A %*% y))
}

REMLlogL <- function(ADcholRinv,
                     ADcholGinv,
                     ADcholC,
                     phi,
                     psi,
                     theta,
                     yPy) {
  # calculute logLikelihood...
  logdetR <- -logdet(ADcholRinv, phi)
  logdetG <- ifelse(!is.null(ADcholGinv), -logdet(ADcholGinv, psi), 0)
  logdetC <- logdet(ADcholC, theta)

  # See e.g. Smith 1995..
  logL <- -0.5 * (logdetR + logdetG + logdetC + yPy)
  return(logL)
}

calcSumSquares <- function(lRinv,
                           Q,
                           r,
                           a,
                           Nvarcomp) {
  SSr <- sapply(X = lRinv, FUN= function(X) {
    quadForm(r, X)
  })
  if (Nvarcomp > 0) {
    SSa <- sapply(X = Q, FUN = function(X) {
      quadForm(a, X)
    })
    SS_all <- c(SSr, SSa)
  } else {
    SS_all <- SSr
  }
  return(SS_all)
}

#' @importFrom stats update
solveMME <- function(cholC,
                     listC,
                     lWtRinvY,
                     phi,
                     theta) {
  C <- linearSum(theta = theta, matrixList = listC)
  cholC <- update(cholC, C)
  WtRinvy <- as.vector(linearSum(theta = phi, matrixList = lWtRinvY))
  a <- spam::backsolve.spam(cholC, spam::forwardsolve.spam(cholC, WtRinvy))
  return(a)
}

sparseMixedModels <- function(y,
                              X,
                              Z,
                              lGinv,
                              lRinv,
                              maxiter = 100,
                              eps = 1.0e-6,
                              display = FALSE,
                              monitor = FALSE) {
  Ntot <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  Nres <- length(lRinv)
  Nvarcomp <- length(lGinv)
  NvarcompTot <- Nres + Nvarcomp
  dimMME <- p + q

  W <- spam::as.spam(cbind(X,Z))
  Wt <- t(W)

  lWtRinvW <- lapply(X = lRinv, FUN = function(x) { Wt %*% x %*% W})
  lWtRinvY <- lapply(X = lRinv, FUN = function(x) { Wt %*% (x %*% y)})
  # extend lGinv with extra zero's for fixed effect...
  lQ <- lapply(X = lGinv, FUN = function(x) {
    zero <- spam::spam(0, ncol = p, nrow = p)
    return(spam::bdiag.spam(zero, x))
  })
  listC <- c(lWtRinvW, lQ)

  # remove some extra zero's....
  lWtRinvW <- lapply(X = lWtRinvW, FUN = spam::cleanup)
  lQ <- lapply(X = lQ, FUN = spam::cleanup)
  lGinv <- lapply(X = lGinv, FUN = spam::cleanup)
  listC <- lapply(X = listC, FUN = spam::cleanup)

  phi <- rep(1.0, Nres)
  psi <- rep(1.0, Nvarcomp)
  theta <- c(psi, phi)
  #cholRinv = chol(Rinv)
  #if (Nvarcomp > 0.0) {cholGinv = chol(linearSum(psi, lGinv))}
  cholC <- chol(linearSum(theta = theta, matrixList = listC))

  if (display) {
    C <- linearSum(theta, listC)
    spam::display(C)
    L <- t(spam::as.spam(cholC))
    spam::display(L)
  }

  # make ADchol for Rinv, Ginv and C:
  ADcholRinv <- ADchol(lRinv)
  if (Nvarcomp > 0) {
    ADcholGinv <- ADchol(lGinv)
  } else {
    ADcholGinv <- NULL
  }
  ADcholC <- ADchol(listC)

  logLprev <- Inf
  if (monitor) {
    cat("iter logLik\n")
  }

  for (it in 1:maxiter) {
    if (Nvarcomp > 0) {
      phi <- theta[c(1:length(phi))]
      psi <- theta[-c(1:length(phi))]
    } else {
      phi <- theta
      psi <- NULL
    }

    # calculate effective dimensions....
    EDmax_phi <- phi * dlogdet(ADcholRinv, phi)
    if (!is.null(ADcholGinv)) {
      EDmax_psi <- psi * dlogdet(ADcholGinv, psi)
    } else {
      EDmax_psi <- NULL
    }
    EDmax <- c(EDmax_phi, EDmax_psi)
    ED <- EDmax - theta * dlogdet(ADcholC, theta)

    # solve mixed model equations and calculate residuals...
    a <- solveMME(cholC = cholC, listC = listC, lWtRinvY = lWtRinvY,
                  phi = phi, theta = theta)
    r <- y - W %*% a

    SS_all <- calcSumSquares(lRinv = lRinv, Q = lQ, r = r, a = a,
                             Nvarcomp = Nvarcomp)
    Rinv <- linearSum(phi, lRinv)
    # here we use Johnson and Thompson 1995:
    yPy <- quadForm(x = y, A = Rinv, y = r)
    # yPy2 = sum(theta*SS_all) (should give same results?)

    logL <- REMLlogL(ADcholRinv = ADcholRinv, ADcholGinv = ADcholGinv,
                     ADcholC = ADcholC, phi = phi, psi = psi, theta = theta,
                     yPy = yPy)

    if (monitor) {
      cat(sprintf("%4d %8.4f\n", it, logL))
    }
    if (abs(logLprev - logL) < eps) {
      break
    }

    theta <- ED / (SS_all + 1.0e-20)
    logLprev <- logL

  }
  names(phi) <- names(lRinv)
  names(psi) <- names(lGinv)
  EDnames <- c(names(lRinv), names(lGinv))
  L <- list(logL = logL, sigma2e = 1.0 / phi, tau2e = 1.0 / psi, ED = ED,
            EDmax = EDmax, EDnames = EDnames, a = a, yhat = W %*% a)
  return(L)
}

#' @importFrom stats model.frame terms model.matrix contrasts
#'
#' @export
LMMsolve <- function(fixed,
                     random = NULL,
                     randomMatrices = NULL,
                     lGinverse = NULL,
                     data,
                     residualterm = NULL,
                     eps = 1.0e-6,
                     monitor = FALSE,
                     display = FALSE,
                     maxiter = 250) {
  ## make random part:
  if (!is.null(random)) {
    mf <- model.frame(random, data, drop.unused.levels = TRUE, na.action = NULL)
    mt <- terms(mf)
    f.terms <- all.vars(mt)[attr(mt,"dataClasses") == "factor"]
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

  if (!is.null(randomMatrices)) {
    ndx <- unlist(randomMatrices)
    dim2.r <- sapply(X = randomMatrices, FUN = length)
    term2.labels.r <- names(randomMatrices)
    Z2 <- data[, ndx]
  } else {
    dim2.r <- NULL
    term2.labels.r <- NULL
    Z2 <- NULL
  }

  if (!(is.null(random) & is.null(randomMatrices))) {
    if (is.null(random)) {
      Z <- Z2
    } else if (is.null(randomMatrices)) {
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

  if (!is.null(residualterm)) {
    lRinv <- makeRlist(df = data, column = residualterm)
  } else {
    lRinv <- list()
    n <- nrow(data)
    lRinv[[1]] <- spam::diag.spam(1, n)
    names(lRinv) = "residual"
  }
  y <- mf[, 1]
  obj <- sparseMixedModels(y, X, Z, lGinv, lRinv, eps, monitor,
                           display = display, maxiter = maxiter)
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



