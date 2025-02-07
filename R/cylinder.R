#' @describeIn spl1D 2-dimensional splines
#'
#' @export
cylinder <- function(x1,
                     x1lim = range(x1),
                  x2,
                  nseg,
                  pord = 2,
                  degree = 3) {
  scaleX <- FALSE
  x2lim <- c(0, 1)
  scaleFactor <- rep(1, 2)

  ## Checks.
  if (!is.numeric(pord) || length(pord) > 1 || pord != 2) {
    stop("pord should be equal to 2 \n")
  }
  if (!is.numeric(degree) || length(degree) > 1 || degree < 1 ||
      degree != round(degree)) {
    stop("degree should be a positive integer.\n")
  }
  if (pord > degree) {
    stop("pord should be less or equal to degree.\n")
  }
  if (!is.numeric(nseg) || length(nseg) != 2 || any(nseg < 1) ||
      any(nseg != round(nseg))) {
    stop("nseg should be a vector of length two of positive integers.\n")
  }
  ## Save names of the x-variables so they can be used later on in predictions.
  x1Name <- deparse(substitute(x1))
  x2Name <- deparse(substitute(x2))
  xNames <- c(x1Name, x2Name)
  missVars <- xNames[!sapply(X = xNames, FUN = exists,
                             where = parent.frame(), inherits = FALSE)]
  if (length(missVars) > 0) {
    stop("The following variables in the spline part of the model ",
         "are not in the data:\n", paste0(missVars, collapse = ", "), "\n",
         call. = FALSE)
  }
  knots <- vector(mode = "list", length = 2)
  knots[[1]] <- PsplinesKnots(x1lim[1], x1lim[2], degree = degree, nseg = nseg[1])
  knots[[2]] <- PsplinesKnots(x2lim[1], x2lim[2], degree = degree, nseg = nseg[2])
  B1 <-  Bsplines(knots[[1]], x1)
  B2 <- cBsplines(knots[[2]], x2)
  q <- c(ncol(B1), ncol(B2))
  B12 <- RowKronecker(B1, B2)

  G1 <- constructG(knots[[1]], scaleX, pord)
  # define null - space
  z2 <- c(0:(nseg[2]-1))/nseg[2]
  u2_x <- sin(2*pi*z2)
  u2_y <- cos(2*pi*z2)
  G2 <- cbind(u2_x, u2_y)

  # check this
  G <- G1[,2] %x% G2
  X1 <- B1 %*% G1
  X2 <- B2 %*% G2
  X12 <- B12 %*% G

  X <- cbind(X1, X2, X12)
  X <- X[,-1]

  ## nominal effective dimension.
  EDnom <- rep(ncol(B12) - ncol(X), 2)

  D1 <- spam::diff.spam(spam::diag.spam(q[1]), diff = pord)
  cD2 <- cDiff(q[2])
  P1 <- spam::crossprod.spam(D1) %x% diag(q[2])
  P2 <- diag(q[1]) %x% spam::crossprod.spam(cD2)

  # add additional constraints, see Boer (2023).
  B1_star <- Bsplines(knots[[1]], x = x1lim)
  cB2_star <- cBsplines(knots[[2]], x = c(0.00, 0.25))
  B_star <- B1_star %x% cB2_star

  P1 <- P1 + crossprod(B_star)
  P2 <- P2 + crossprod(B_star)
  lGinv <- list(P1 = P1, P2 = P2)
  names(lGinv) <- paste0("s(", xNames, ")")
  dim.f <- ncol(X)
  term.labels.f <- paste0("fix(", paste(xNames, collapse = ", "), ")")
  dim.r <- ncol(B12)
  term.labels.r <-  paste0("s(", paste(xNames, collapse = ", "), ")")
  xList <- setNames(list(x1, x2), xNames)
  return(list(X = X, Z = B12, lGinv = lGinv, knots = knots,
              dim.f = dim.f, dim.r = dim.r, term.labels.f = term.labels.f,
              term.labels.r = term.labels.r, x = xList, pord = pord,
              degree = degree, scaleX = scaleX, EDnom = EDnom,
              scaleFactor = scaleFactor))
}

