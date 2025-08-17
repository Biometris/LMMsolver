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
  SSr2 <- sapply(X = lWtRinvY, FUN = function(x) { sum(a * x) })
  SSr3 <- sapply(X = lWtRinvW, FUN = function(x) { quadForm(a, x)})
  SSr <- SSr1 - 2 * SSr2 + SSr3
  if (Nvarcomp > 0) {
    SSa <- sapply(X = lQ, FUN = function(x) {
      quadForm(a, x)
    })
    SS_all <- c(SSa, SSr)
  } else {
    SS_all <- SSr
  }
  # make sure sum of squares are all positive.
  SS_all <- pmax(SS_all, 1e-12)
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
                              tolerance = 1e-6,
                              trace = FALSE,
                              theta = NULL,
                              fixedTheta = NULL,
                              grpTheta = NULL) {
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
      xW@entries <- xW@entries * rep(d, times = diff(xW@rowpointers))
    } else {
      xW <- x %*% W
    }
    MatrixProduct(t(W), xW)
  })
  lWtRinvY <- lapply(X = lRinv, FUN = function(x) {
    spam::crossprod.spam(W, x %*% y) })
  lYtRinvY <- lapply(X = lRinv, FUN= function(x) {
    quadForm(y, x) })

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
  lRinv <- lapply(X = lRinv, FUN = spam::cleanup)
  lC <- lapply(X = lC, FUN = spam::cleanup)
  if (is.null(theta)) {
    theta <- rep(1, Nvarcomp + Nres)
    theta_restr <- theta
  } else {
    theta_restr <- theta
  }
  if (is.null(fixedTheta)) {
    ## Fix a penalty theta, if value becomes high.
    if (!is.null(grpTheta)) {
      fixedTheta <- rep(FALSE, length = max(grpTheta))
    } else {
      fixedTheta <- rep(FALSE, length = NvarcompTot)
    }
  }
  if (!is.null(grpTheta)) {
    nGrp <- max(grpTheta)
    if (length(fixedTheta) != nGrp) {
      stop("problem with number of groups defined in grpTheta argument", call. = FALSE)
    }
    conM <- spam::spam(x = 0, nrow = NvarcompTot, ncol = nGrp)
    for (i in seq_len(NvarcompTot)) {
      conM[i, grpTheta[i]] <- 1
    }
    fixedThetaRes <- fixedTheta
  } else {
    conM <- spam::diag(1, NvarcompTot)
    fixedThetaRes <- fixedTheta
  }

  ## Check the stucture of Rinv, don't allow for overlapping penalties
  M <- sapply(lRinv, FUN = function(x) {
    abs(spam::diag.spam(x)) > getOption("spam.eps")})
  rSums <- rowSums(M)
  #if (!isTRUE(all.equal(rSums, rep(1, nrow(M))))) {
  if (max(rSums) > 1) {
    stop("Overlapping penalties for residual part.\n", call. = FALSE)
  }
  ## Calculate number of elements per group for residuals.
  nR <- colSums(M)
  EDmax_phi <- nR
  Rinv <- Reduce("+", lRinv)
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
  traceDf <- NULL
  for (it in seq_len(maxit)) {
    if (Nvarcomp > 0) {
      psi <- theta[seq_len(Nvarcomp)]
      phi <- theta[-(seq_len(Nvarcomp))]
    } else {
      psi <- NULL
      phi <- theta
    }
    ## calculate logdet for Rinv.
    logdetR <- -sum(nR * log(phi)) - logdetRinvConstant

    ## calculated logdet and dlogdet for Ginv and C.
    ## Ginv, if exists
    if (!is.null(ADcholGinv)) {
      dlogdetGinv <- dlogdet(ADcholGinv, psi)
      logdetG <- -attr(dlogdetGinv, which = "logdet")
    } else {
      logdetG <- 0
    }
    ## update the expressions including Rinv.
    YtRinvY <- sum(phi * unlist(lYtRinvY))
    WtRinvY <- as.vector(linearSum(theta = phi, matrixList = lWtRinvY))

    ## matrix C.
    dlogdetC <- dlogdet(ADcholC, theta, WtRinvY)
    logdetC <- attr(dlogdetC, which = "logdet")
    a <- attr(dlogdetC, which = "x.coef")

    ## calculate effective dimensions.
    if (!is.null(ADcholGinv)) {
      EDmax_psi <- psi * dlogdetGinv
    } else {
      EDmax_psi <- NULL
    }
    EDmax <- c(EDmax_psi, EDmax_phi)
    ED <- EDmax - theta * dlogdetC

    ## to make sure ED is always positive.
    ED <- pmax(ED, .Machine$double.eps)

    ## calculate Sum of Squares.
    SS_all <- calcSumSquares(lYtRinvY, lWtRinvY, lWtRinvW, lQ, a, Nvarcomp)

    ## Johnson and Thompson 1995, see below eq [A5]: Py = R^{-1} r.
    yPy <- YtRinvY - sum(a * WtRinvY)

    ## calculate REMLlogL, see e.g. Smith 1995.
    logL <- -0.5 * (logdetR + logdetG + logdetC + yPy)
    ## Add to trace.
    traceDf <- rbind(traceDf, c(iter = it, logLik = logL, ED))

    if (trace) {
      cat(sprintf("%4d %8.4f\n", it, logL))
    }
    if (abs(logLprev - logL) < tolerance) {
      break
    }
    ED_restr <- as.vector(ED %*% conM)
    SS_all_restr <- as.vector(SS_all %*% conM)
    theta_restr <- ifelse(fixedThetaRes, theta_restr, ED_restr / SS_all_restr)
    theta <- as.vector(conM %*% theta_restr)

    fixedTheta <- theta > 1.0e6 | fixedTheta
    fixedThetaRes <- theta_restr > 1.0e6 | fixedThetaRes

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
  traceDf <- data.frame(traceDf)
  colnames(traceDf)[-(1:2)] <- paste("ED", EDnames)
  L <- list(logL = logL, sigma2e = 1 / phi, tau2e = 1 / psi, ED = ED,
            theta = theta, EDnames = EDnames, a = a, yhat = yhat,
            residuals = r, nIter = it, C = C, cholC = cholC,
            trace = traceDf)
  return(L)
}

