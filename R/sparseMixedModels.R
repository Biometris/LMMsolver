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
  sum(x * (A %*% y))
}

REMLlogL <- function(ADcholRinv,
                     ADcholGinv,
                     ADcholC,
                     phi,
                     psi,
                     theta,
                     yPy) {

  # calculate logLikelihood...
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
    SS_all <- c(SSa, SSr)
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


#' @keywords internal
sparseMixedModels <- function(y,
                              X,
                              Z,
                              lGinv,
                              lRinv,
                              maxit = 100,
                              tolerance = 1.0e-6,
                              trace = FALSE) {
  Ntot <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  Nres <- length(lRinv)
  Nvarcomp <- length(lGinv)
  NvarcompTot <- Nres + Nvarcomp
  dimMME <- p + q

  W <- spam::as.spam(cbind(X, Z))
  Wt <- t(W)

  lWtRinvW <- lapply(X = lRinv, FUN = function(x) { Wt %*% x %*% W})
  lWtRinvY <- lapply(X = lRinv, FUN = function(x) { Wt %*% (x %*% y)})
  # extend lGinv with extra zero's for fixed effect...
  lQ <- lapply(X = lGinv, FUN = function(x) {
    zero <- spam::spam(0, ncol = p, nrow = p)
    return(spam::bdiag.spam(zero, x))
  })
  listC <- c(lQ, lWtRinvW)

  # remove some extra zero's....
  lWtRinvW <- lapply(X = lWtRinvW, FUN = spam::cleanup)
  lQ <- lapply(X = lQ, FUN = spam::cleanup)
  lGinv <- lapply(X = lGinv, FUN = spam::cleanup)
  listC <- lapply(X = listC, FUN = spam::cleanup)

  psi <- rep(1.0, Nvarcomp)
  phi <- rep(1.0, Nres)
  theta <- c(psi, phi)

  C <- linearSum(theta = theta, matrixList = listC)
  opt <- summary(C)
  cholC <- chol(C, memory = list(nnzR = 8 * opt$nnz,
                                 nnzcolindices = 4 * opt$nnz))

  # make ADchol for Rinv, Ginv and C:
  ADcholRinv <- ADchol(lRinv)
  if (Nvarcomp > 0) {
    ADcholGinv <- ADchol(lGinv)
  } else {
    ADcholGinv <- NULL
  }
  ADcholC <- ADchol(listC)

  logLprev <- Inf
  if (trace) {
    cat("iter logLik\n")
  }

  for (it in 1:maxit) {
    if (Nvarcomp > 0) {
      psi <- theta[c(1:length(psi))]
      phi <- theta[-c(1:length(psi))]
    } else {
      psi <- NULL
      phi <- theta
    }

    # calculate effective dimensions....
    EDmax_phi <- phi * dlogdet(ADcholRinv, phi)
    if (!is.null(ADcholGinv)) {
      EDmax_psi <- psi * dlogdet(ADcholGinv, psi)
    } else {
      EDmax_psi <- NULL
    }
    EDmax <- c(EDmax_psi, EDmax_phi)
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

    if (trace) {
      cat(sprintf("%4d %8.4f\n", it, logL))
    }
    if (abs(logLprev - logL) < tolerance) {
      break
    }

    theta <- ED / (SS_all + 1.0e-20)
    logLprev <- logL

  }
  if (it == maxit) {
    warning("No convergence after ", maxit, " iterations \n", call. = FALSE)
  }
  C <- linearSum(theta = theta, matrixList = listC)
  cholC <- update(cholC, C)

  names(phi) <- names(lRinv)
  names(psi) <- names(lGinv)
  EDnames <- c(names(lGinv), names(lRinv))
  yhat <- W %*% a

  L <- list(logL = logL, sigma2e = 1 / phi, tau2e = 1 / psi, ED = ED,
            theta = theta,
            EDmax = EDmax, EDnames = EDnames, a = a, yhat = yhat,
            residuals = y - yhat, nIter = it, C = C)
  return(L)
}




