# helper to convert from matrix to vec form
mat2vec <- function(C, Mask) {
  len <- nrow(Mask)
  theta <- numeric(len)
  for (k in seq_len(len)) {
    i <- Mask[k, 1]
    j <- Mask[k, 2]
    if (i==j) {
      theta[k] <- C[i, j]
    } else {
      theta[k] <- 2*C[i, j]
    }
  }
  theta
}

# helper to convert from vec to matrix form
vec2mat <- function(theta, Mask)
{
  d <- max(Mask)
  len <- nrow(Mask)
  C <- matrix(0, d, d)
  for (k in seq_len(len)) {
    i <- Mask[k, 1]
    j <- Mask[k, 2]
    if (i==j) {
      C[i, i] <- theta[k]
    } else {
      C[i, j] <- C[j,i] <- 0.5*theta[k]
    }
  }
  C
}

# forward function eta -> xi
eta_to_xi <- function(eta, Mask) {
  d <- max(Mask)
  len <- nrow(Mask)
  M <- matrix(data=0, nrow=d, ncol=d)
  for (k in seq_len(len)) {
    i <- Mask[k, 1]
    j <- Mask[k, 2]
    if (i==j) {
      M[i,j] <- exp(0.5*eta[k])
    } else {
      M[i,j] <- sinh(eta[k])*M[i,i]
    }
  }
  M
}

# forward function xi -> theta
xi_to_theta <- function(xi, Mask) {
    thetaMat <- tcrossprod(xi)
    mat2vec(thetaMat, Mask)
}

# backward function theta -> xi
theta_to_xi <- function(M) {
  t(chol(M))
}

# backward function xi -> eta
xi_to_eta <- function(xi, Mask)
{
  d <- max(Mask)
  len <- nrow(Mask)
  eta <- numeric(length=len)
  for (k in seq_len(len)) {
    i <- Mask[k, 1]
    j <- Mask[k, 2]
    if (i==j) {
      eta[k] <- 2*log(xi[i,i])
    } else {
      eta[k] <- asinh(xi[i,j]/xi[i,i])
    }
  }
  eta
}

# rotate in theta - space
rotate_theta <- function(M, xi) {
  crossprod(xi, M %*% xi)
}

gradient <- function(eta, Mask, y, U, ADC, ADP, ADRinv, lC) {

  UtY <- t(U) %*% y

  xi <- eta_to_xi(eta, Mask)

  theta <- xi_to_theta(xi, Mask)
  K <- length(theta) - 1

  #
  # begin calculations in theta (precision space)
  #
  phi <- theta[K+1]
  psi <- theta[-(K+1)]
  dlogdetC <- dlogdet(ADC, theta, phi*UtY)
  dlogdetP <- dlogdet(ADP, psi)
  dlogdetRinv <- dlogdet(ADRinv, phi)
  a <- attr(dlogdetC,"x.coef")
  r <- y - U %*% a
  ss <- rep(NA, K+1)
  dlogdet <- rep(NA, K+1)
  for (k in 1:K) {
    ss[k] <- quadForm(a, lC[[k]])
    dlogdet[k] <- dlogdetP[k] - dlogdetC[k]
  }
  ss[K+1] <- sum(r^2)
  dlogdet[K+1] <- dlogdetRinv[1] - dlogdetC[K+1]   #

  ED <- theta*dlogdet
  Tmat <- vec2mat(dlogdet, Mask)

  Smat <- vec2mat(ss, Mask)

  T_rotated <- rotate_theta(Tmat, xi)
  S_rotated <- rotate_theta(Smat, xi)
  #
  # end calculations in theta (precision) space.
  #

  xi_T <- theta_to_xi(T_rotated)
  xi_S <- theta_to_xi(S_rotated)

  eta_T <- xi_to_eta(xi_T, Mask)
  eta_S <- xi_to_eta(xi_S, Mask)

  grad <- eta_T - eta_S
  L <- list (grad = grad,
             theta = theta,
             a = a,
             r = r,
             phi = phi,
             ED = ED,
             dlogdetC = dlogdetC,
             dlogdetP = dlogdetP,
             dlogdetRinv = dlogdetRinv)
  L
}

HarvilleODE <- function(y, X, Z, lGinv, lRinv, Mask, alpha, maxiter, thr)
{
  p <- ncol(X)
  U <- spam::as.spam(cbind(X, Z))
  UtU <- crossprod(U)
  UtY <- t(U) %*% y
  K <- length(lGinv)
  lC <- list()
  for (k in 1:K) {
    lC[[k]] <- spam::bdiag.spam(spam::diag.spam(0,p), lGinv[[k]])
  }
  lC[[K+1]] <- UtU

  # construct ADchol from a list of semi-positive precision matrices
  ADC <- ADchol(lC)
  ADP <- ADchol(lGinv)
  ADRinv <- ADchol(lRinv)
  eta <- rep(0, K+1)
  logPrev <- Inf
  logL <- 0
  trace <- NULL
  for (it in seq_len(maxiter)) {
    z <- gradient(eta, Mask, y, U, ADC, ADP, ADRinv, lC)

    grad <- z$grad
    a <- z$a
    ED <- z$ED
    r <- z$r
    phi <- z$phi

    eta <- eta + alpha*grad

    yPy <- sum(y*r)*phi
    logdetRinv <- attr(z$dlogdetRinv, which ='logdet')
    logdetPinv <- attr(z$dlogdetP, which='logdet')
    logdetC <- attr(z$dlogdetC, which='logdet')
    logL <- 0.5 * (logdetRinv + logdetPinv - logdetC - yPy)

    names(eta) <- paste0("eta", 1:(K+1))
    names(ED) <- paste0("ED", 1:(K+1))
    trace <- rbind(trace,
                   data.frame(it=it, t(eta), t(ED) , logL = logL))

    if (abs(logL-logPrev) < thr) break

    logPrev <- logL
  }
  L <- list(a = a,
            theta = z$theta,
            niter = it,
            ED = ED,
            logL = logL,
            trace = trace)
  L
}



