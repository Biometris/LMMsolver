
h <- function(eta) {
  v <- exp(eta)
  v/(1+sum(v))
}

Jacobian <- function(eta) {
  nc <- length(eta)
  v <- exp(eta)
  g <- 1+sum(v)
  J <- matrix(data=0, nrow=nc, ncol=nc)
  for (r in seq_len(nc)) {
    for (j in seq_len(nc)) {
      J[r,j] <- -(v[r]/g)*(v[j]/g)
      if (r==j) {
        J[r,r] <- J[r,r] + v[r]/g
      }
    }
  }
  return(J)
}

calcSigma <- function(pi) {
  nc <- length(pi)
  Sigma <- matrix(data=0, nrow=nc, ncol=nc)
  for (r in seq_len(nc)) {
    for (j in seq_len(nc)) {
      if (r==j) {
        Sigma[r, r] <- pi[r]*(1-pi[r])
      } else {
        Sigma[r, j] <- -pi[r]*pi[j]
      }
    }
  }
  return(Sigma)
}

#' Defines the multinomial model
#'
#' @export
multinomial <- function() {
  family <- "multinomial"
  structure(list(family = family), class = "family")
}

#' Solver for multinomial distribution, at the moment only for 1D splines
#'
#' @export
Multinomialsolve <- function(fixed, dat, family = multinomial()) {

  family <- multinomial()
  if (family$family != "multinomial") {
    stop("Not multinomial()!")
  }

  mf <- model.frame(fixed, dat, drop.unused.levels = TRUE, na.action = NULL)
  Y <- model.response(mf, type = "any")
  colNamesResponse <- colnames(Y)
  # nCategories (minus one)
  nCat <- length(colNamesResponse) - 1
  sY <- rowSums(Y)
  Y <- Y/sY
  range(rowSums(Y))
  dat_fr <- data.frame(x=x,Y)
  Y <- Y[,-(nCat+1)] # remove last column

  y <- as.vector(t(Y))
  x <- dat$x

  knots <- PsplinesKnots(0, 1, degree = 3, nseg = 50)
  B <- Bsplines(knots, x)
  X <- cbind.spam(1, x) %x% diag(nCat)
  Z <- B %x% diag(nCat)
  q <- ncol(B)
  scaleFactor <- calcScaleFactor(list(knots), pord=2)
  lGinv <- constructGinvSplines(q, list(knots), pord=2, scaleFactor)
  lGinv[[1]] <- lGinv[[1]] %x% diag(nCat)

  # remark: We can add some prior information to initial
  # values based on the scores:
  eta <- rep(0, n*nCat)

  # Initialize the sparse block structure for W and Dinv
  W <- diag.spam(1, n) %x% spam(x=1, nrow=nCat, ncol=nCat)
  Dinv <- diag.spam(1, n) %x% spam(x=1, nrow=nCat, ncol=nCat)

  theta <- c(1,1)
  fixedTheta <- c(FALSE, TRUE)
  for (it in 1:100) {
    Eta <- matrix(data=eta, ncol=nCat, nrow=n, byrow=TRUE)
    Pi <- t(apply(X=Eta, MARGIN=1, FUN=h))
    pi <- as.vector(t(Pi))

    D_list <- lapply(1:n, FUN =
                       function(r) {
                         return(Jacobian(Eta[r,]))
                       } )
    Dinv_list <- lapply(1:n, FUN= function(r) { solve(D_list[[r]])})

    W_list <- lapply(1:n, FUN =
                       function(r) {
                         Sigma <- calcSigma(Pi[r,])
                         D_list[[r]] %*% solve(Sigma) %*% t(D_list[[r]])
                       } )

    # Store in the defined sparse format for W and Dinv:
    Dinv@entries <- unlist(Dinv_list)
    W@entries <- unlist(W_list)

    z <- eta + Dinv %*% (y - pi)
    lRinv <- list(W)
    obj <- sparseMixedModels(z, X = X, Z = Z,
                             lGinv = lGinv, lRinv = lRinv, trace=TRUE,
                             fixedTheta = fixedTheta, theta = theta)
    eta_old  <- eta
    eta <- obj$yhat
    theta <- obj$theta
    tol <- sum((eta - eta_old)^2)/sum(eta^2)
    cat("iter ", it, "     ", tol, "    ", obj$logL, "   ", obj$ED, "\n")
    if (tol < 1.0e-10) {
      cat("convergence after", it, "iterations\n")
      break;
    }
  }
  return(list(eta = eta,
               colNamesResponse = colNamesResponse,
              knots = knots,
              obj = obj))
}
