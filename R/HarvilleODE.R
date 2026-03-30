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

HarvilleODE <- function(y, X, Z, lGinv, lRinv, Mask, alpha, maxiter, thr=1.0e-6, trace=FALSE)
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
  df.trace <- NULL
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
    df.trace <- rbind(df.trace,
                   data.frame(it=it, t(eta), t(ED) , logL = logL))
    if (trace) {
      cat("it", it, logL, "\n")
    }

    if (abs(logL-logPrev) < thr) break

    logPrev <- logL
  }
  C <- linearSum(z$theta, lC)
  L <- list(a = a,
            theta = z$theta,
            niter = it,
            ED = ED,
            logL = logL,
            C = C,
            trace = df.trace)
  L
}


# help function
# this only works for cubic B-splines,
# it replaces the builtin-constraints Boer2023
# with hard constraints.
sparse_NP_mat <- function(q) {
  # N_P: q x (q-2)
  N <- matrix(0, nrow = q, ncol = q - 2)

  # Fill identity for "free" coefficients (a2 ... a_{q-1})
  for (i in 1:(q-2)) {
    N[i+1, i] <- 1
  }

  # Left boundary: a1 = -4 a2 - a3
  N[1, 1] <- -4
  if (q > 3) {
    N[1, 2] <- -1
  }

  # Right boundary: a_q = -a_{q-2} - 4 a_{q-1}
  N[q, q-2] <- -4
  if (q > 3) {
    N[q, q-3] <- -1
  }

  # Return sparse
  spam::as.spam(N)
}


# simple test function to model genotypic specific curves, including Random Regression
fit_subject_specific <-function(data, nseg, maxiter=250,thr=1.0e-6,trace=FALSE) {

  deg <- 3
  pord <- 2
  # assume for the moment that data has columns y, geno as factor, and time:
  y <- data$y
  geno <- data$geno
  time <- data$time
  n <- nrow(data)
  nGeno <- nlevels(geno)
  X <- cbind(1, time)
  # construction of knots for P-splines
  knots <- PsplinesKnots(xmin = min(time), max(time), degree=deg, nseg=nseg)
  B <- Bsplines(knots, time)
  q <- ncol(B)
  Np <- sparse_NP_mat(q)

  Z_spl <- B %*% Np

  # Construct Penalty, with built-in boundary constraints:
  D <- spam::diff.spam(spam::diag.spam(q), diff=2)
  DtD <- crossprod(D)
  B_ref <- Bsplines(knots, x=c(min(time),max(time)))
  range(B_ref %*% Np)
  P <- t(Np) %*% DtD %*% Np

  D_G <- t(spam::diff.spam(spam::diag.spam(nGeno), diff=1))
  Z_geno <- model.matrix(~geno-1, data = data)
  Z_geno_D_G <- Z_geno %*% D_G

  Z_con_geno <- RowKronecker(Z_geno_D_G, Z_spl)

  ran_regr <- RowKronecker(Z_geno, as.spam(time))

  Z <- cbind(Z_geno, ran_regr, Z_spl, Z_con_geno)

  P_con_geno <- crossprod(D_G) %x% P

  lM <- list()
  lM[[1]] <- tcrossprod(c(1, 0))
  lM[[2]] <- tcrossprod(c(0, 1))
  lM[[3]] <- spam::spam(x=c(0, 0.5, 0.5, 0), nrow = 2, ncol = 2)
  K <- length(lM)
  lM_ext <- list()
  for (k in 1:K) {
    lM_ext[[k]] <- spam::as.spam(lM[[k]] %x% spam::diag.spam(1, nGeno))
  }

  lGinv <- list()
  lGinv[[1]] <- spam::bdiag.spam(lM_ext[[1]], spam::diag.spam(0, q-2), spam::diag.spam(0,(nGeno-1)*(q-2)))
  lGinv[[2]] <- spam::bdiag.spam(lM_ext[[2]], spam::diag.spam(0, q-2), spam::diag.spam(0,(nGeno-1)*(q-2)))
  lGinv[[3]] <- spam::bdiag.spam(lM_ext[[3]], spam::diag.spam(0, q-2), spam::diag.spam(0,(nGeno-1)*(q-2)))
  lGinv[[4]] <- spam::bdiag.spam(spam::diag.spam(0, 2*nGeno), P, spam::diag.spam(0,(nGeno-1)*(q-2)))
  lGinv[[5]] <- spam::bdiag.spam(spam::diag.spam(0, 2*nGeno), spam::diag.spam(0, q-2), P_con_geno)

  lRinv <- list(Rinv = spam::diag.spam(n))

  Mask <- matrix(data=c(1,1 , 2,2, 2,1, 3,3, 4,4, 5,5), nrow = 6, ncol = 2, byrow=TRUE)
  obj <- HarvilleODE(y, X, Z, lGinv, lRinv, Mask,
                                 alpha=1.0, maxiter = maxiter, thr = thr,
                                 trace=trace)

  # calculate vcov
  P <- linearSum(obj$theta[1:3], lM)
  Sigma <- solve(P)
  Sigma

  L <- list(model = obj,
            knots = knots,
            Sigma = Sigma,
            Np    = Np,
            q     = q,
            D_G   = D_G)
  L
}

predict_subject_specific <- function(object, nGrid) {
  obj <- object$model
  knots <- object$knots
  Np <- object$Np
  q <- object$q
  D_G <- object$D_G
  xmin <- attr(knots, which="xmin")
  xmax <- attr(knots, which="xmax")

  time_grid <- seq(xmin, xmax, length=nGrid)
  Bg <- LMMsolver:::Bsplines(knots, time_grid)

  p <- 2
  beta <- obj$a[1:p]
  ndx_mu_ran <- p + geno
  ndx_beta_ran <- nGeno + p + geno
  mu_ran <- obj$a[ndx_mu_ran]
  beta_ran <- obj$a[ndx_beta_ran]
  non_lin_con <- obj$a[-c(1:(2*nGeno+2+q-2))]
  coef_mean <- obj$a[(2*nGeno+2+1):(2*nGeno+2+q-2)]

  grid <- expand.grid(
    time = time_grid,
    geno = factor(1:nGeno)
  )

  Bg <- Bsplines(knots, grid$time)
  BgNp <- Bg %*% Np

  grid$mu_ran   <- mu_ran[as.numeric(grid$geno)]
  grid$beta_ran <- beta_ran[as.numeric(grid$geno)]

  grid$f_lin <- beta[1] + beta[2]*grid$time +
    grid$mu_ran + grid$beta_ran * grid$time

  grid$f_mean <- as.vector(BgNp %*% coef_mean)
  grid$f_mean_tot <- grid$f_mean + beta[1] + beta[2]*grid$time

  coef_dev <- matrix(non_lin_con, nrow = (nGeno-1), byrow = TRUE)

  # Map genotypes to contrasts
  Dg_mat <- as.matrix(D_G)

  Zg_dev <- Dg_mat[as.numeric(grid$geno), ]

  grid$f_dev <- rowSums((Zg_dev %*% coef_dev) * as.matrix(BgNp))

  grid$fit <- grid$f_lin + grid$f_mean + grid$f_dev

  grid
}





