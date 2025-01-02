
# generalized logit, see Fahrmeir et al.
glogit <- function(eta) {
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

