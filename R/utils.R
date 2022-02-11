#' Norm of a vector
#'
#' @param x a numerical vector.
#'
#' @return norm of vector x
#'
#' @keywords internal
normVec <- function(x) {
  return(sqrt(sum(x ^ 2)))
}

#' Combine two lists of Ginv matrices
#'
#' output is a combined list, of dimension (d1+d2) x (d1+d2),
#' with matrices in lGinv1 extended to [A1 0]
#'                                     [0  0]
#' and              lGinv2 extended to [0  0]
#'                                     [0 A2]
#' @return a list of sparse matrices of dimension (d1+d2) x (d1+d2)
#'
#' @param lGinv1 a list of sparse matrices.
#' @param lGinv2 a list of sparse matrices.
#'
#' @keywords internal
expandGinv <- function(lGinv1, lGinv2) {
  if (is.null(lGinv1)) return(lGinv2)
  if (is.null(lGinv2)) return(lGinv1)

  d1 <- nrow(lGinv1[[1]])
  d2 <- nrow(lGinv2[[1]])

  zero1 <- spam::diag.spam(0, d2)
  zero2 <- spam::diag.spam(0, d1)

  lGinv1a <- lapply(X = lGinv1, FUN = function(x) {
    spam::bdiag.spam(x, zero1)
  })
  lGinv2a <- lapply(X = lGinv2, FUN = function(x) {
    spam::bdiag.spam(zero2, x)
  })
  lGinv <- c(lGinv1a, lGinv2a)
  return(lGinv)
}

#' Construct P-splines penalty matrix D'D
#'
#' @param q integer with dimensions.
#' @param pord order of the penalty.
#'
#' @return qxq matrix D'D of class spam
#'
#' @keywords internal
constructPenalty <-function(q,
                            pord) {
  D <- spam::diff.spam(spam::diag.spam(q), diff = pord)
  DtD <- spam::crossprod.spam(D)
  return(DtD)
}

#' Construct fixed part of the spline model
#'
#' @param B matrix with B-spline basis.
#' @param x vector with values for x.
#' @param scaleX logical. If scaleX is FALSE, the original x is used. If scaleX
#' is TRUE, scaling is used, based on the B-splines basis. For details see
#' the code.
#' @param pord order of the penalty, values 1 or 2.
#'
#' @return a matrix X
#'
#' @keywords internal
constructX <- function(B,
                       x,
                       scaleX,
                       pord) {
  if (pord == 2) {
    if (scaleX) {
      ## calculate the linear/fixed parts.
      U_null <- cbind(1, scale(seq_len(ncol(B))))
      U_null <- apply(U_null, MARGIN = 2, function(x) (x / normVec(x)))
      X <- B %*% U_null
    } else {
      X <- cbind(1, x)
    }
  } else {  # pord = 1
    X <- matrix(data = 1, nrow = length(x), ncol = 1)
  }
  X
}

#' Construct constraint matrix
#'
#' @param q dimension of the B-spline basis used.
#' @param pord order of the penalty matrix (pord=1 or 2).
#'
#' @return a q x q matrix of type spam
#'
#' @keywords internal
constructCCt <- function(q,
                         pord) {
  C <- spam::spam(x = 0, nrow = q, ncol = pord)
  C[1, 1] <- C[q, pord] <- 1
  CCt <- C %*% t(C)
  return(CCt)
}


#' Construct list of Ginv matrices of splines
#'
#' @param q vector, with dimension of the B-spline basis used.
#' @param pord order of the penalty matrix (1 or 2).
#'
#' @return a list of symmetric matrices of length of vector q.
#'
#' @keywords internal
constructGinvSplines <- function(q,
                                 pord) {
  ## dimension
  d <- length(q)
  lCCt <- lapply(X = q, FUN = constructCCt, pord = pord)
  CCt <- Reduce("kronecker", lCCt)
  lGinv <- list()
  for (i in seq_len(d)) {
    L <- list()
    for (j in seq_len(d)) {
      if (i == j) {
        L[[j]] <- constructPenalty(q[j], pord = pord)
      } else {
        L[[j]] <- spam::diag.spam(q[j])
      }
    }
    lGinv[[i]] <- Reduce("kronecker", L) + CCt
  }
  return(lGinv)
}


#' Remove the intercept from a design matrix
#'
#' @param X design matrix.
#'
#' @return a matrix if \code{X} has more than one column, otherwise return NULL
#'
#' @keywords internal
removeIntercept <- function(X) {
  if (ncol(X) == 1) {
    X <- NULL
  } else {
    X <- X[, -1,  drop = FALSE]
  }
  return(X)
}

#' Calculate the Nominal Effective dimension
#'
#' @keywords internal
calcNomEffDim <- function(X,
                          Z,
                          dim.r) {
  if (is.null(Z)) return(NULL)
  p <- ncol(X)
  EDnom <- vector(length = length(dim.r))
  e <- cumsum(dim.r)
  s <- e - dim.r + 1
  # for each variance component in Z:
  for (i in 1:length(dim.r)) {
    ndx <- s[i]:e[i]
    Zi <- Z[, ndx]
    # if number of columns is high, use upper bound:
    if (dim.r[i] > 100) {
      colSum <- colSums(Zi)
      if (all(colSum==colSum[1]))
        EDnom[i] <- dim.r[i] - 1
      else
        EDnom[i] <- dim.r[i]
    } else {
      XZ <- cbind(X, Zi)
      r <- qr(XZ)$rank
      EDnom[i] <- r - p
    }
  }
  return(EDnom)
}

#' Check variables in formula
#'
#' Check that all variables in a formula are present in the data.
#'
#' @importFrom utils hasName
#'
#' @keywords internal
checkFormVars <- function(formula,
                          data) {
  if (!is.null(formula)) {
    ## Get variables in formula.
    formVars <- all.vars(formula)
    missVars <- formVars[!hasName(data, formVars)]
    if (length(missVars) > 0) {
      ## Get the name of formula as passed to the function.
      formName <- deparse(substitute(formula))
      stop("The following variables in the ", formName, " part of the model ",
           "are not in the data:\n", paste0(missVars, collapse = ", "), "\n",
           call. = FALSE)
    }
  }
}

#' @noRd
grp <- function(x) {
  NULL
}

#' Check group part
#'
#' Check that all variables in group are also specified as grp() in random part
#' of the model and vice versa.
#'
#' @importFrom utils hasName
#'
#' @keywords internal
checkGroup <- function(random,
                       group) {
  ## Only check if at least one of random and group is not NULL.
  if (!is.null(random) || !is.null(group)) {
    if (is.null(random)) {
      grpVars <- NULL
    } else {
      ## Extract variables specified in grp() functions in random part.
      randTerms <- terms(random, specials = "grp")
      grpPos <- attr(randTerms, which = "specials")$grp
      grpTerms <- sapply(X = grpPos, FUN = function(pos) {randTerms[pos]})
      grpVars <- sapply(X = grpTerms, FUN = all.vars)
      ## Check for variables in grp missing in group.
      if (is.null(group)) {
        grpMiss <- grpVars
      } else {
        grpMiss <- grpVars[sapply(X = grpVars, FUN = function(grpVar) {
          !hasName(x = group, name = grpVar)
        })]
      }
      if (length(grpMiss) > 0) {
        stop("The following variables are specified in grp in the random part ",
             "of the model but not present in group:\n",
             paste0(grpMiss, collapse = ", "), "\n", call. = FALSE)
      }
      ## Remove grp terms from random part of the model.
      if (length(grpPos) > 0) {
        random <- randTerms[-grpPos]
      }
    }
    ## Remove variables from group missing in grp.
    if (any(names(group) %in% grpVars)) {
      group <- group[names(group) %in% grpVars]
    } else {
      group <- NULL
    }
  }
  return(list(random = random, group = group))
}

#' Helper function for naming coefficients.
#'
#' @keywords internal
nameCoefs <- function(coefs,
                      termLabels,
                      s,
                      e,
                      data,
                      group = NULL,
                      type = "fixed") {
  coefRes <- vector(mode = "list", length = length(termLabels))
  for (i in seq_along(coefRes)) {
    coefI <- coefs[s[i]:e[i]]
    labI <- termLabels[i]
    ## Split labI in interaction terms (if any).
    labISplit <- strsplit(x = labI, split = ":")[[1]]
    if (length(labISplit) == 1) {
      ## No interactions.
      if (labI == "(Intercept)") {
        names(coefI) <- "(Intercept)"
      } else if (startsWith(x = labI, prefix = "lin(") ||
                 startsWith(x = labI, prefix = "s(")) {
        ## Spline terms are just named 1...n.
        names(coefI) <- paste0(labI, "_", seq_along(coefI))
      } else if (!is.null(group) && labI %in% names(group)) {
        ## For group combine group name and column name.
        names(coefI) <- paste(labI, colnames(data)[group[[labI]]], sep = "_")
      } else if (is.factor(data[[labI]])) {
        ## For fixed terms an extra 0 for the reference level has to be added.
        if (type == "fixed") {
          coefI <- c(0, coefI)
        }
        names(coefI) <- paste(labI, levels(data[[labI]]), sep = "_")
      } else {
        ## Numerical variable. Name equal to label.
        names(coefI) <- labI
      }
    } else {
      ## Interaction term.
      if (all(sapply(X = labISplit, FUN = function(labTerm) {
        is.factor(data[[labTerm]])
      }))) {
        namesTerm <- lapply(X = labISplit, FUN = function(labTerm) {
          paste(labTerm, levels(data[[labTerm]]), sep = "_")
        })
        ## For fixed terms an extra 0 for the reference level has to be added.
        if (type == "fixed") {
          coefI <- c(0, coefI)
        }
        names(coefI) <- levels(interaction(namesTerm, sep = ":"))
      } else {
        ## Numerical variable. Name equal to label.
        names(coefI) <- labI
      }
    }
    coefRes[[i]] <- coefI
  }
  return(coefRes)
}


