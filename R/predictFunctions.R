#' Standard errors for predictions
#'
#' Calculates the standard errors for predictions \eqn{D \hat{u}},
#' see Welham et al. 2004 and Gilmour et al. 2004 for details.
#'
#' The prediction error variance is given by \eqn{D C^{-1} D'},
#' where \eqn{C} is the mixed model coefficient matrix, and \eqn{D} defines
#' linear combinations of fixed and random effects.
#' The standard errors are given by the the square root of
#' the diagonal. To calculate the standard errors in an efficient way we use that
#'
#' \deqn{\frac{\partial log|C + \xi_i d_i d_i'|}{\partial \xi_i} |_{\xi_i=0}
#'  = trace(C^{-1} d_i d_i') =
#' trace(d_i' C^{-1} d_i) = d_i' C^{-1} d_i, }
#' where \eqn{d_i} is row \eqn{i} of matrix \eqn{D}. The values of
#' \eqn{d_i' C^{-1} d_i} can be calculated more efficient, avoiding the
#' calculation of the inverse of \eqn{C}, by using Automated Differentiation
#' of the Choleksy algorithm, see section 2.3 in Smith (1995) for details.
#'
#' @param C a symmetric matrix of class spam
#' @param D a matrix of class spam
#'
#' @returns a vector with standard errors for predictions \eqn{D \hat{u}}.
#'
#' @references
#' Welham, S., Cullis, B., Gogel, B., Gilmour, A., & Thompson, R. (2004).
#' Prediction in linear mixed models.
#' Australian & New Zealand Journal of Statistics, 46(3), 325-347.
#'
#' @references
#' Smith, S. P. (1995). Differentiation of the Cholesky algorithm.
#' Journal of Computational and Graphical Statistics, 4(2), 134-147.
#'
#' @references
#' Gilmour, A., Cullis, B., Welham, S., Gogel, B., & Thompson, R. (2004).
#' An efficient computing strategy for prediction in mixed linear models.
#' Computational statistics & data analysis, 44(4), 571-586.
#'
#' @keywords internal
calcStandardErrors <- function(C,
                               D) {
  tD <- spam::t.spam(D)
  DtD <- tD %*% D

  ## It adds extra zeros ("fill-ins") to matrix C, needed
  ## to calculate the Partial Derivatives of Cholesky.
  iDtD <- DtD
  iC <- C
  iC@entries <- rep(2, length(iC@entries))
  iDtD@entries <- rep(1, length(iDtD@entries))
  Cwf <- iC + iDtD
  inC <- which(Cwf@entries==2 | Cwf@entries == 3)
  Cwf@entries <- rep(0, length(Cwf@entries))
  Cwf@entries[inC] <- C@entries

  cholC <- spam::chol.spam(Cwf)
  p <- cholC@pivot
  tDp <- tD[p, ]
  x <- diagXCinvXt(cholC, tDp)
  se <- sqrt(x)
  return(se)
}

auxFun <- function(object, i, nr, nc, Names=NULL)
{
  # MB: More efficient to use t(Dg), working row-wise:
  Dg <- spam::spam(x = 0, nrow = nr, ncol = nc)
  ndx <- object$ndxCoefficients[[i]]
  len_ndx <- length(ndx)
  rm <- which(ndx==0)
  if (length(rm)>0) {
    ndx <- ndx[-which(ndx==0)]
  }
  if (!is.null(Names)) {
    for (j in seq_along(ndx)) {
      Dg[,ndx[j]] <- 1*(Names==names(ndx)[j])
    }
  } else { # averaging
    Dg[,ndx] <- 1/len_ndx
  }
  Dg
}


#' Test function for predict, for the moment internal
#'
#' @keywords internal
predictTest <- function(object, classify) {
  # count what kind of interactions, if any:
  cat("Predict function: \n")
  s <- classify
  s2 <- gsub(":","", s)
  numOcc <- nchar(s) - nchar(s2)
  numOcc

  if (numOcc == 1) {
    Interaction <- TRUE
  } else {
    Interaction <- FALSE
  }

  if (Interaction) {
    # names of marginals
    marg_names <- unlist(strsplit(classify, split=":"))
  } else {
    marg_names <- NULL
  }

  term.labels <- c(object$term.labels.f, object$term.labels.r)
  isFixed <- term.labels %in% object$term.labels.f
  dim <- object$dim
  e <- cumsum(dim)
  s <- e - dim + 1
  df.terms <- data.frame(label = term.labels, s = s, e = e, fixed = isFixed)

  ndxPred <- object$ndxCoefficients[[classify]]

  namesI <- names(ndxPred)

  nr <- length(ndxPred)
  nc <- ncol(object$C)

  ndx_intercept <- which(df.terms$label == '(Intercept)')
  ndx_classify <- which(df.terms$label == classify)
  ndx_marginals <- which(df.terms$label %in% marg_names)
  ndx_averaging <- which(df.terms$fixed)
  ndx_averaging <- setdiff(ndx_averaging, ndx_intercept)
  ndx_averaging <- setdiff(ndx_averaging, ndx_classify)
  ndx_averaging <- setdiff(ndx_averaging, ndx_marginals)
  ndx_averaging

  cat("  Intercept: ", df.terms$label[ndx_intercept],"\n")
  cat("  Classify:  ", df.terms$label[ndx_classify], "\n")
  cat("  Marginals: ", df.terms$label[ndx_marginals], "\n")
  cat("  Averaging: ", df.terms$label[ndx_averaging], "\n")

  # intercept:
  Dg <- spam::spam(x = 0, nrow = nr, ncol = nc)
  Dg[, 1] <- 1
  # ndx for classify:
  Dg <- Dg + auxFun(object, ndx_classify, nr, nc, namesI)
  for (i in ndx_marginals) {
    label <- df.terms$label[i]          # label of current term
    k <- which(label==marg_names)
    namesLab <- sapply(namesI, FUN = function(x) { strsplit(x,split=":")[[1]][k] })
    Dg <- Dg + auxFun(object, i, nr, nc, namesLab)
  }
  for (i in ndx_averaging) {
    Dg <- Dg + auxFun(object, i, nr, nc, NULL)
  }
  Dg <- spam::cleanup(Dg)

  ypred <- as.vector(Dg %*% object$coefMME)
  ypredse <- calcStandardErrors(object$C, Dg)

  if (Interaction) {
    marg1_labels <- sapply(namesI, FUN = function(x) { strsplit(x,split=":")[[1]][1]})
    marg2_labels <- sapply(namesI, FUN = function(x) { strsplit(x,split=":")[[1]][2]})

    marg1 <- sapply(marg1_labels, FUN=function(x) { strsplit(x,"_")[[1]][2]})
    marg2 <- sapply(marg2_labels, FUN=function(x) { strsplit(x,"_")[[1]][2]})
    pred <- data.frame(marg1, marg2, ypred, ypredse)
    colnames(pred) <- c(marg_names, "prediction", "se")
    # order on first marginal
    ord <- order(pred[[marg_names[1]]])
    pred <- pred[ord, ]
  } else {
    term <- sapply(namesI, FUN=function(x) { strsplit(x,"_")[[1]][2]})
    pred <- data.frame(term, ypred, ypredse)
    colnames(pred) <- c(classify, "prediction", "se")
  }

  rownames(pred) <- NULL
  pred
}



