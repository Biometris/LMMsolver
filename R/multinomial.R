
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

extend_coef <- function(ndx, respVar) {

  nCat <- length(respVar)-1
  Cat <- respVar[-c(nCat+1)]

  for (i in seq_along(ndx)) {
    tType <- attr(ndx[[i]], which="termType")
    x <- as.numeric(ndx[[i]])
    y <- as.vector(sapply(x,FUN=function(x) {
      if (x==0) return(rep(0,nCat))
      nCat*(x-1) + c(1:nCat)}))
    namesCoef <- names(ndx[[i]])
    name_ext <- rep(namesCoef, each=nCat)
    cat_ext <- rep(Cat, times=length(x))
    name_ext <- paste0(name_ext, "_", cat_ext)
    names(y) <- name_ext
    ndx[[i]] <- y
    attr(ndx[[i]],"termType") <- tType
  }
  ndx
}



#' Defines the multinomial model
#'
#' @export
multinomial <- function() {
  family <- "multinomial"
  structure(list(family = family), class = "family")
}

