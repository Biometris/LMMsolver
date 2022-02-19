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
  WtRinvy <- as.vector(linearSum(theta = phi, matrixList = lWtRinvY))
  a <- spam::backsolve.spam(cholC, spam::forwardsolve.spam(cholC, WtRinvy))
  return(a)
}

#' @keywords internal
calcEffDim <- function(ADcholGinv,
                       ADcholRinv,
                       ADcholC,
                       phi,
                       psi,
                       theta) {
  dlogdetRinv <- dlogdet(ADcholRinv, phi)
  logdetR <- -attr(dlogdetRinv, which = "logdet")

  # Ginv, if exists:
  if (!is.null(ADcholGinv)) {
    dlogdetGinv <- dlogdet(ADcholGinv, psi)
    logdetG <- -attr(dlogdetGinv, "logdet")
  } else {
    logdetG <- 0
  }
  ## matrix C.
  dlogdetC <- dlogdet(ADcholC, theta)
  logdetC <- attr(dlogdetC, which = "logdet")
  ## calculate effective dimensions.
  EDmax_phi <- phi * dlogdetRinv
  if (!is.null(ADcholGinv)) {
    EDmax_psi <- psi * dlogdetGinv
  } else {
    EDmax_psi <- NULL
  }
  EDmax <- c(EDmax_psi, EDmax_phi)
  ED <- EDmax - theta * dlogdetC
  attributes(ED) <- NULL
  attr(ED, "logdetG") <- logdetG
  attr(ED, "logdetC") <- logdetC
  attr(ED, "logdetR") <- logdetR
  return(ED)
}

#' @keywords internal
REMLlogL <- function(ED,
                     yPy) {
  logdetC <- attr(ED, "logdetC")
  logdetG <- attr(ED, "logdetG")
  logdetR <- attr(ED, "logdetR")
  # REML-loglikelihood (without constant), see e.g. Smith 1995.
  logL <- -0.5 * (logdetR + logdetG + logdetC + yPy)
  return(logL)
}


#' @keywords internal
sparseMixedModels <- function(y,
                              X,
                              Z,
                              lGinv,
                              lRinv,
                              maxit = 100,
                              tolerance = 1.0e-6,
                              trace = FALSE,
                              theta = NULL) {
  Ntot <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  Nres <- length(lRinv)
  Nvarcomp <- length(lGinv)
  NvarcompTot <- Nres + Nvarcomp
  dimMME <- p + q
  W <- spam::as.spam(cbind(X, Z))
  lWtRinvW <- lapply(X = lRinv, FUN = function(x) {
    spam::crossprod.spam(W, x %*% W) })
  lWtRinvY <- lapply(X = lRinv, FUN = function(x) {
    spam::crossprod.spam(W, x %*% y) })
  ## Extend lGinv with extra zeros for fixed effect.
  lQ <- lapply(X = lGinv, FUN = function(x) {
    zero <- spam::spam(0, ncol = p, nrow = p)
    return(spam::bdiag.spam(zero, x))
  })
  listC <- c(lQ, lWtRinvW)
  ## Remove some extra zeros.
  lWtRinvW <- lapply(X = lWtRinvW, FUN = spam::cleanup)
  lQ <- lapply(X = lQ, FUN = spam::cleanup)
  lGinv <- lapply(X = lGinv, FUN = spam::cleanup)
  listC <- lapply(X = listC, FUN = spam::cleanup)
  if (is.null(theta)) {
    theta <- rep(1, Nvarcomp + Nres)
  }
  if (Nvarcomp > 0) {
    psi <- theta[1:Nvarcomp]
    phi <- theta[-(1:Nvarcomp)]
  } else {
    psi <- NULL
    phi <- theta
  }
  C <- linearSum(theta = theta, matrixList = listC)
  opt <- summary(C)
  cholC <- chol(C, memory = list(nnzR = 8 * opt$nnz,
                                 nnzcolindices = 4 * opt$nnz))
  ## Make ADchol for Rinv, Ginv and C:
  ADcholRinv <- ADchol(lRinv)
  if (Nvarcomp > 0) {
    ADcholGinv <- ADchol(lGinv)
  } else {
    ADcholGinv <- NULL
  }
  ADcholC <- ADchol(listC)
  ## Initialize values for loop.
  logLprev <- Inf
  ## Fix a penalty theta, if value becomes high.
  fixedTheta <- rep(FALSE, length = NvarcompTot)
  if (trace) {
    cat("iter logLik\n")
  }
  for (it in 1:maxit) {
    if (Nvarcomp > 0) {
      psi <- theta[1:length(psi)]
      phi <- theta[-(1:length(psi))]
    } else {
      psi <- NULL
      phi <- theta
    }
    ## calculate the effective dimension (plus logdet as attributes).
    ED <- calcEffDim(ADcholGinv, ADcholRinv, ADcholC, phi, psi, theta)

    ## update the cholesky with new parameters theat
    C <- linearSum(theta = theta, matrixList = listC)
    cholC <- update(cholC, C)

    ## solve mixed model equations.
    a <- solveMME(cholC = cholC, listC = listC, lWtRinvY = lWtRinvY,
                  phi = phi, theta = theta)
    ## calculate the residuals.
    r <- y - W %*% a
    SS_all <- calcSumSquares(lRinv = lRinv, Q = lQ, r = r, a = a,
                             Nvarcomp = Nvarcomp)
    Rinv <- linearSum(phi, lRinv)
    ## Johnson and Thompson 1995, see below eq [A5]: Py = R^{-1} r
    yPy <- quadForm(x = y, A = Rinv, y = r)
    ## calculate REMLlogL, ED has logdet as attributes:
    logL <- REMLlogL(ED, yPy)
    if (trace) {
      cat(sprintf("%4d %8.4f\n", it, logL))
    }
    if (abs(logLprev - logL) < tolerance) {
      break
    }
    ## Update the penalties theta that are not fixed.
    theta <- ifelse(fixedTheta, theta, ED / SS_all)
    ## Set elements of theta fixed if penalty > 1.0e6.
    fixedTheta <- theta > 1.0e6
    logLprev <- logL
  }
  if (it == maxit) {
    warning("No convergence after ", maxit, " iterations \n", call. = FALSE)
  }
  names(phi) <- names(lRinv)
  names(psi) <- names(lGinv)
  EDnames <- c(names(lGinv), names(lRinv))
  L <- list(logL = logL, sigma2e = 1 / phi, tau2e = 1 / psi, ED = ED,
            theta = theta, EDnames = EDnames, a = a, yhat = y - r,
            residuals = r, nIter = it, C = C, cholC = cholC)
  return(L)
}

