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

# build from the restrictions and the penalties an indicator matrix
build_restriction_matrix <- function(C_restrict, lP) {
  n_penalty <- length(lP)
  n_constraints <- ncol(C_restrict)
  kappa <- matrix(0, nrow = n_penalty, ncol = n_constraints)
  for (i in seq_len(n_penalty)) {
    for (j in seq_len(n_constraints)) {
      cj <- C_restrict[, j]
      kappa[i, j] <- as.numeric(t(cj) %*% lP[[i]] %*% cj)
    }
  }
  kappa <- 1*(kappa > 1.0e-12)
  kappa
}

PrecisionMatrix_lifting <- function(C_restrict, lP, kappa) {
  n_penalties <- length(lP)
  n_constraints <- ncol(C_restrict)

  lP_lifted <- lP

  for (i in seq_len(n_penalties)) {

    # select only constraints relevant for this penalty
    idx <- which(kappa[i, ] == 1)

    if (length(idx) > 0) {
      C_i <- C_restrict[, idx, drop = FALSE]

      lP_lifted[[i]] <- lP[[i]] + C_i %*% spam::t.spam(C_i)
    }
  }

  lP_lifted
}


## ---- log-det contribution (psi-dependent part) ----
logdet_correction <- function(M, psi)
{
  logdetQ <- 0
  for (j in seq_len(ncol(M))) {
    idx <- which(M[, j] == 1)
    logdetQ <- logdetQ + log(sum(psi[idx]))
  }
  logdetQ
}

## ---- EDs per penalty ----
ED_corrections <- function(M, psi) {
  ED <- numeric(length(psi))
  for (j in seq_len(ncol(M))) {
    idx <- which(M[, j] == 1)
    ED[idx] <- ED[idx] + psi[idx] / sum(psi[idx])
  }
  ED
}

# make the matrix positive definite
# PrecisionMatrix_lifting <- function(C_restrict, lP) {
#   n_penalties <- length(lP)
#   lP_lifted <- lP
#   for (i in seq_len(n_penalties)) {
#     C_i <- C_restrict[, i, drop = FALSE]
#     lP_lifted[[i]] <- lP[[i]] + C_i %*% spam::t.spam(C_i)
#   }
#   lP_lifted
# }

compute_alpha_active <- function(x, y, active, lower = -7, upper = 7) {
  alpha_max <- Inf

  for (i in seq_along(x)) {
    if (!active[i]) next  # skip fixed variables

    if (y[i] > 0) {
      alpha_i <- (upper - x[i]) / y[i]
      alpha_max <- min(alpha_max, alpha_i)
    } else if (y[i] < 0) {
      alpha_i <- (lower - x[i]) / y[i]
      alpha_max <- min(alpha_max, alpha_i)
    }
  }

  return(alpha_max)
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
                              grpTheta = NULL,
                              C_restrict = NULL) {
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

  if (!is.null(C_restrict)) {
    kappa <- build_restriction_matrix(C_restrict, lGinv)
    lGinv <- PrecisionMatrix_lifting(C_restrict, lGinv, kappa)
  }

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
      if (!is.null(C_restrict)) {
        logdetG <- logdetG + logdet_correction(kappa, psi)
      }
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
      if (!is.null(C_restrict)) {
        EDmax_psi <- EDmax_psi - ED_corrections(kappa, psi)
      }
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

    # do updat on a log-scale
    eta <- log10(theta_restr)
    grad_eta <- log10(ED_restr/theta_restr) - log10(SS_all_restr)

    # freeze gradients
    grad_eta[fixedThetaRes] <- 0

    alpha_max <- compute_alpha_active(
      eta, grad_eta, !fixedThetaRes,
      lower = -6, upper = 6
    )

    alpha <- min(1, alpha_max)

    eta_new <- eta + alpha * grad_eta

    # detect boundary hits
    tol <- 1e-10
    hit_upper <- eta_new >= (6 - tol)
    hit_lower <- eta_new <= (-6 + tol)

    eta_new[hit_upper] <- 6
    eta_new[hit_lower] <- -6

    # update active set
    fixedThetaRes <- fixedThetaRes | hit_upper | hit_lower

    #if (all(fixedThetaRes)) {
    #  stop("Cannot find solution; all penalties fixed\n")
    #}

    #cat("eta ", eta, "  eta_new", eta_new, " non-active " , fixedThetaRes, "\n")

    # back-transform
    theta_restr <- 10^eta_new
    theta <- as.vector(conM %*% theta_restr)

    #fixedTheta <- theta > 1.0e6 | fixedTheta
    #fixedThetaRes <- theta_restr > 1.0e6 | fixedThetaRes

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

