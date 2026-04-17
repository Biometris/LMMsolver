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
fit_subject_specific <-function(data, nseg, nseg.ad, maxiter=250,thr=1.0e-6,trace=FALSE) {

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

  ran_regr <- RowKronecker(Z_geno, spam::as.spam(time))

  Z <- cbind(Z_geno, ran_regr, Z_spl, Z_con_geno)

  P_con_geno <- crossprod(D_G) %x% P

  if (nseg.ad > 0) {
    # adaptive, for ridge penalty :
    knots.ad <- LMMsolver:::PsplinesKnots(xmin=1,xmax=q,nseg=nseg.ad,degree=1)
    C.ad <- LMMsolver:::Bsplines(knots.ad, x=1:q)
  }

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

  if (nseg.ad > 0) {
    # adaptive part
    dim_C.ad <- ncol(C.ad)
    cat("Dim C.ad ", dim(C.ad), "\n")
    cat("Dim Np   ", dim(Np), "\n")
    for (i in seq_len(dim_C.ad)) {
      P_r <- (t(Np) %*% spam::diag.spam(C.ad[,i]) %*% Np)
      P_ridge <- crossprod(D_G) %x% P_r
      lGinv[[5+i]] <- spam::bdiag.spam(spam::diag.spam(0, 2*nGeno), spam::diag.spam(0, q-2), P_ridge)
    }
  }

  lRinv <- list(Rinv = spam::diag.spam(n))

  Mask <- matrix(data=c(1,1 , 2,2, 2,1, 3,3, 4,4, 5,5), nrow = 6, ncol = 2, byrow=TRUE)
  if (nseg.ad > 0) {
    for (i in seq_len(dim_C.ad)) {
      Mask <- rbind(Mask, c(5+i,5+i))
    }
  }

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
            nGeno = nGeno,
            D_G   = D_G)
  L
}

predict_subject_specific <- function(object, nGrid) {
  nGeno <- object$nGeno
  geno <- 1:nGeno
  obj <- object$model
  knots <- object$knots
  Np <- object$Np
  q <- object$q
  D_G <- object$D_G
  xmin <- attr(knots, which="xmin")
  xmax <- attr(knots, which="xmax")

  time_grid <- seq(xmin, xmax, length=nGrid)
  Bg <- Bsplines(knots, time_grid)

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

  grid$f_lin_dev <- grid$mu_ran + grid$beta_ran * grid$time


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





