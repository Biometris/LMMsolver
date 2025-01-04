# generalized logit, see Fahrmeir et al.
glogit <- function(mu) {
  log(mu/(1-sum(mu)))
}

# generalized inverse logit, see Fahrmeir et al.
inv_glogit <- function(eta) {
  v <- exp(eta)
  v/(1+sum(v))
}

x_logx <- function(x) {
  nonzero <- (x > 1.0e-20)
  z <- x
  z[nonzero] <- x[nonzero]*log(x[nonzero])
  z
}

multinomial.deviance.residuals <- function(y, mu, wt) {
  devy <- x_logx(y)
  devmu <- x_logx(mu)
  dev <- 2 * wt * (devy - devmu)
  return(dev)
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

#' Family Object for Multinomial Model
#'
#' The Multinomial model is not part of the standard family. The implementation
#' is based on Chapter 6 in Fahrmeir et al. (2013).
#'
#' @references Fahrmeir, Ludwig, Thomas Kneib, Stefan Lang, Brian Marx, Regression models.
#' Springer Berlin Heidelberg, 2013.
#'
#' @returns
#' An object of class \code{familyLMMsolver} with the following components:
#' \item{family}{character string with the family name.}
#' \item{linkfun}{the link function.}
#' \item{linkinv}{the inverse of the link function.}
#' \item{dev.resids}{function giving the deviance for each observation as a function of (y, mu, wt)}
#'
#' @export
multinomial <- function() {
  family <- "multinomial"
  link <- glogit
  linkinv <- inv_glogit
  dev.resids <- multinomial.deviance.residuals
  structure(list(family = family,
                 link = link,
                 linkinv = linkinv,
                 dev.resids = dev.resids), class = "familyLMMsolver")
}

