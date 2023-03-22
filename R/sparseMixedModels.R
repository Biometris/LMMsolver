# calcalate the weighted sum of list of matrices.
linearSum  <- function(theta,
                       matrixList) {
  C <- Reduce('+', Map(f = function(x, y) x * y, theta, matrixList))
  return(C)
}

# calculate quadratic form x'Ay (efficient):
quadForm <- function(x,
                     A,
                     y = x) {
  sum(x * (A %*% y))
}

calcSumSquares <- function(lYtRinvY,
                           lWtRinvY,
                           lWtRinvW,
                           lQ,
                           a,
                           Nvarcomp) {
  SSr1 <- unlist(lYtRinvY)
  SSr2 <- sapply(X=lWtRinvY, FUN = function(X) { sum(a*X) })
  SSr3 <- sapply(X=lWtRinvW, FUN = function(X) { quadForm(a, X)})
  SSr <- SSr1 - 2*SSr2 + SSr3
  if (Nvarcomp > 0) {
    SSa <- sapply(X = lQ, FUN = function(X) {
      quadForm(a, X)
    })
    SS_all <- c(SSa, SSr)
  } else {
    SS_all <- SSr
  }
  # make sure sum of squares are all positive.
  SS_all <- sapply(X = SS_all, FUN = function(x) {max(1.0e-12, x)})
  return(SS_all)
}

#' @importFrom stats update
#' @keywords internal
sparseMixedModels <- function(y,
                              X,
                              Z,
                              lGinv,
                              lRinv,
                              maxit = 100,
                              tolerance = 1.0e-6,
                              trace = FALSE,
                              theta = NULL,
                              fixedTheta = NULL) {
  Ntot <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  Nres <- length(lRinv)
  Nvarcomp <- length(lGinv)
  NvarcompTot <- Nres + Nvarcomp
  dimMME <- p + q
  W <- spam::cbind.spam(X, Z)
  lWtRinvW <- lapply(X = lRinv, FUN = function(x) {
    d <- spam::diag.spam(x)
    # if x is diagonal calculate x %*% W in a more efficient way
    if (isTRUE(all.equal(spam::diag.spam(d), x))) {
      xW <- W
      xW@entries <- xW@entries * rep(d, times=diff(xW@rowpointers))
    } else {
      xW <- x %*% W
    }
    spam::crossprod.spam(W, xW)
  })
  lWtRinvY <- lapply(X = lRinv, FUN = function(x) {
    spam::crossprod.spam(W, x %*% y) })
  lYtRinvY <- lapply(X =lRinv, FUN= function(x) {
    quadForm(y,x) })

  ## Extend lGinv with extra zeros for fixed effect.
  lQ <- lapply(X = lGinv, FUN = function(x) {
    zero <- spam::spam(0, ncol = p, nrow = p)
    return(spam::bdiag.spam(zero, x))
  })
  lC <- c(lQ, lWtRinvW)
  ## Remove some extra zeros.
  lWtRinvW <- lapply(X = lWtRinvW, FUN = spam::cleanup)
  lQ <- lapply(X = lQ, FUN = spam::cleanup)
  lGinv <- lapply(X = lGinv, FUN = spam::cleanup)
  lC <- lapply(X = lC, FUN = spam::cleanup)
  if (is.null(theta)) {
    theta <- rep(1, Nvarcomp + Nres)
  }
  if (is.null(fixedTheta)) {
    ## Fix a penalty theta, if value becomes high.
    fixedTheta <- rep(FALSE, length = NvarcompTot)
  }
  if (Nvarcomp > 0) {
    psi <- theta[1:Nvarcomp]
    phi <- theta[-(1:Nvarcomp)]
  } else {
    psi <- NULL
    phi <- theta
  }

  ## Check the stucture of Rinv, don't allow for overlapping
  ## penalties
  M <- sapply(lRinv,FUN=function(x) {
    as.integer(abs(spam::diag.spam(x))>getOption("spam.eps"))})
  rSums <- rowSums(M)
  #if (!isTRUE(all.equal(rSums, rep(1, nrow(M))))) {
  if (max(rSums) > 1) {
    stop("overlapping penalties for residual part") }

  # calculate number of elements per group for residuals
  nR <- colSums(M)
  EDmax_phi <- nR
  Rinv <- Reduce("+",lRinv)
  logdetRinvConstant <- as.numeric(spam::determinant.spam(Rinv)$modulus)

  ## Make ADchol for Ginv and C:
  if (Nvarcomp > 0) {
    ADcholGinv <- ADchol(lGinv)
  } else {
    ADcholGinv <- NULL
  }
  ADcholC <- ADchol(lC)
  ## Initialize values for loop.
  logLprev <- Inf
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

    ## calculate logdet for Rinv
    logdetR <- -sum(nR*log(phi)) - logdetRinvConstant

    ## calculated logdet and dlogdet for Ginv and C
    ## Ginv, if exists
    if (!is.null(ADcholGinv)) {
      dlogdetGinv <- dlogdet(ADcholGinv, psi)
      logdetG <- -attr(dlogdetGinv, which = "logdet")
    } else {
      logdetG <- 0
    }
    ## update the expressions including Rinv
    YtRinvY <- sum(phi * unlist(lYtRinvY))
    WtRinvY <- as.vector(linearSum(theta = phi, matrixList = lWtRinvY))

    ## matrix C
    dlogdetC <- dlogdet(ADcholC, theta, WtRinvY)
    logdetC <- attr(dlogdetC, which = "logdet")
    a <- attr(dlogdetC,which="x.coef")

    ## calculate effective dimensions.
    if (!is.null(ADcholGinv)) {
      EDmax_psi <- psi * dlogdetGinv
    } else {
      EDmax_psi <- NULL
    }
    EDmax <- c(EDmax_psi, EDmax_phi)
    ED <- EDmax - theta * dlogdetC

    ## to make sure ED is always positive
    ED <- pmax(ED, .Machine$double.eps)

    ## calculate Sum of Squares
    SS_all <- calcSumSquares(lYtRinvY, lWtRinvY, lWtRinvW, lQ, a, Nvarcomp)

    ## Johnson and Thompson 1995, see below eq [A5]: Py = R^{-1} r
    yPy <- YtRinvY - sum(a*WtRinvY)

    ## calculate REMLlogL, see e.g. Smith 1995
    logL <- -0.5 * (logdetR + logdetG + logdetC + yPy)

    if (trace) {
      cat(sprintf("%4d %8.4f\n", it, logL))
    }
    if (abs(logLprev - logL) < tolerance) {
      break
    }
    ## Update the penalties theta that are not fixed.
    theta <- ifelse(fixedTheta, theta, ED / SS_all)
    ## Set elements of theta fixed if penalty > 1.0e6.
    fixedTheta <- (theta > 1.0e6) | (fixedTheta == TRUE)
    logLprev <- logL
  }
  if (it == maxit) {
    warning("No convergence after ", maxit, " iterations \n", call. = FALSE)
  }

  ## MB: not really needed, just to keep consistent with previous versions.
  C <- linearSum(theta = theta, matrixList = lC)
  opt <- summary(C)
  cholC <- chol(C, memory = list(nnzR = 8 * opt$nnz,
                                 nnzcolindices = 4 * opt$nnz))

  ## calculate yhat and residuals
  yhat <- W %*% a
  r <- y - yhat
  names(phi) <- names(lRinv)
  names(psi) <- names(lGinv)
  EDnames <- c(names(lGinv), names(lRinv))
  L <- list(logL = logL, sigma2e = 1 / phi, tau2e = 1 / psi, ED = ED,
            theta = theta, EDnames = EDnames, a = a, yhat = yhat,
            residuals = r, nIter = it, C = C, cholC = cholC)
  return(L)
}

