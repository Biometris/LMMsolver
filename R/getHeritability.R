## Internal: spectral generalized heritability for a random term
## Returns full spectral decomposition with attribute "h2_G"
.spectralHeritability <- function(obj, geno.term, tol = 1e-8) {

  ## --- indices ---------------------------------------------------------------
  p <- ncol(obj$X)

  ndx_fix  <- seq_len(p)
  ndx_rand <- unlist(obj$ndxCoefficients)
  ndx_rand <- ndx_rand[ndx_rand > 0]

  ndx_all  <- unique(c(ndx_fix, ndx_rand))
  ndx_g    <- as.numeric(obj$ndxCoefficients[[geno.term]])
  ndx_g    <- ndx_g[ndx_g > 0]
  ndx_nuis <- setdiff(ndx_all, ndx_g)

  ## --- genetic covariance ----------------------------------------------------
  i <- which(obj$term.labels.r == geno.term)
  sigma2_g <- 1 / obj$theta[i]

  Ginv_all <- obj$lGinv[[geno.term]]
  ndx_loc  <- ndx_g - p          # ONLY needed here
  Ginv <- Ginv_all[ndx_loc, ndx_loc]
  G <- sigma2_g * solve(Ginv)

  ## --- absorb nuisance effects (Johnson & Thompson) ---------------------------
  C <- obj$C
  C_gg <- C[ndx_g, ndx_g]

  if (length(ndx_nuis) > 0) {
    C_nuis <- C[ndx_nuis, ndx_nuis]
    C_ng   <- C[ndx_nuis, ndx_g]
    C_gn   <- C[ndx_g, ndx_nuis]
    C_gg_abs <- C_gg - C_gn %*% solve(C_nuis, C_ng)
  } else {
    C_gg_abs <- C_gg
  }

  ## --- spectral decomposition (genetic space only) ----------------------------
  m <- ncol(G)

  one <- rep(1, m)
  dG <- solve(G, one)

  PG <- diag(m) - tcrossprod(dG) / drop(t(dG) %*% dG)
  eig_PG <- eigen(PG, symmetric = TRUE)
  U <- eig_PG$vectors[, eig_PG$values > tol, drop = FALSE]

  GZtPZG <- G - solve(C_gg_abs)

  A <- t(U) %*% GZtPZG %*% U
  B <- t(U) %*% G %*% U

  eigB <- eigen(B, symmetric = TRUE)
  B_inv_sqrt <- eigB$vectors %*%
    diag(1 / sqrt(eigB$values)) %*%
    t(eigB$vectors)

  S <- B_inv_sqrt %*% A %*% B_inv_sqrt
  eigS <- eigen(S, symmetric = TRUE)

  lambda <- eigS$values
  Q <- eigS$vectors

  rho <- diag(t(Q) %*% B %*% Q)
  w <- rho / sum(rho)

  out <- data.frame(
    component = seq_along(lambda),
    lambda    = lambda,
    w         = w,
    h2_comp   = w * lambda
  )

  attr(out, "h2_G") <- sum(out$h2_comp)
  out
}

#' Generalized heritability of a random term
#'
#' Computes the generalized heritability of a random-effect term from a fitted
#' linear mixed model. By default, a single scalar heritability value is returned.
#' Optionally, a spectral decomposition is provided that reveals how genetic
#' signal is distributed across estimable genetic directions.
#'
#' @param obj An object of class \code{"LMMsolve"}.
#' @param geno.term A character string giving the name of the genetic random-effect
#'   term.
#' @param type Character string specifying the output:
#'   \code{"scalar"} (default) returns a single numeric heritability value;
#'   \code{"spectral"} returns a data frame with the spectral decomposition.
#' @param tol Numerical tolerance used to determine the estimable genetic space.
#'
#' @return
#' If \code{type = "scalar"}, a numeric value giving the generalized heritability.
#' If \code{type = "spectral"}, a data frame with columns:
#' \describe{
#'   \item{component}{Index of the spectral component}
#'   \item{lambda}{Canonical heritability for the component}
#'   \item{w}{Weight of the component (genetic capacity)}
#'   \item{h2_comp}{Contribution of the component to total heritability}
#' }
#' In the spectral case, the scalar generalized heritability is also available as
#' the attribute \code{"h2_G"}.
#'
#' @details
#' Generalized heritability is defined as the proportion of estimable genetic signal
#' retained by the design relative to the available genetic capacity, accounting
#' for the genetic covariance structure.
#'
#' For independent genotypes, this definition reduces to classical generalized
#' heritability measures based on effective dimension (Cullis, Oakey,
#' Rodríguez-Álvarez). When genotypes are correlated, the spectral decomposition
#' reveals anisotropy in information retention across genetic directions.
#'
#' @export
getHeritability <- function(obj,
                            geno.term,
                            type = c("scalar", "spectral"),
                            tol = 1e-8) {

  ## --- existing checks (unchanged, for tinytest compatibility) ----------------
  if (!inherits(obj, "LMMsolve")) {
    stop("obj must be an object of class 'LMMsolve'")
  }

  if (!is.character(geno.term) || length(geno.term) != 1L) {
    stop("geno.term must be a single character string")
  }

  EDdf <- effDim(obj)

  if (!(geno.term %in% EDdf$Term)) {
    stop(paste(geno.term, "not defined in the model"))
  }

  penalty <- subset(EDdf, Term == geno.term)$Penalty
  if (penalty < 1.0e-15) {
    stop(paste(geno.term, "should be in the random term"))
  }

  ## --- new spectral implementation -------------------------------------------
  type <- match.arg(type)

  spec <- .spectralHeritability(obj, geno.term, tol = tol)

  if (type == "scalar") {
    return(attr(spec, "h2_G"))
  } else {
    return(spec)
  }
}

