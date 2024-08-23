#' Norm of a vector
#'
#' @param x a numerical vector.
#'
#' @return norm of vector x
#'
#' @noRd
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
#' @noRd
#' @keywords internal
expandGinv <- function(lGinv1,
                       lGinv2) {
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

#' Calculate scale factor for precision matrices.
#'
#' To make the penalty matrix more stable if there are many knots,
#' a scaled version is used, \eqn{(1/dx)^(2pord-1) D'D}.
#'
#' @noRd
#' @keywords internal
calcScaleFactor <- function(knots,
                            pord) {
  dx <- sapply(X = knots, FUN = attr, which = "dx")
  sc <- (1 / dx)^(2 * pord - 1)
  sc <- ifelse(sc < 1e-10, 1e-10, sc)
  return(sc)
}

#' Construct the penalty matrix
#'
#' Construct a scaled version of the P-splines penalty matrix, see details.
#'
#' The P-spline penalty matrix has the form \eqn{D'D}, where \eqn{D} is the
#' `pord` order difference matrix.
#'
#' @param q integer with dimensions.
#' @param pord order of the penalty.
#' @param dx distance between the knots.
#'
#' @return qxq penalty matrix of class spam
#'
#' @noRd
#' @keywords internal
constructPenalty <- function(q,
                             pord) {
  D <- spam::diff.spam(spam::diag.spam(q), diff = pord)
  DtD <- spam::crossprod.spam(D)
  return(DtD)
}

#' Consruct Greville points,
#'
#' The Greville points \eqn{\xi_{1,j,d}} are defined by
#' \deqn{\xi_{i,j,p} = \frac{x_{j+1} + \ldots + x_{j+d}}{d}}
#' see Tom Lyche et al.
#' @noRd
#' @keywords internal
GrevillePoints <- function(knots) {
  d <- attr(knots,"degree")
  q <- length(knots) - (d + 1)
  sapply(X = c(1:q), FUN = function(j) { sum(knots[(j+1):(j+d)]) / d })
}

#' Construct G, such that BG = X
#'
#' @noRd
#' @keywords internal
constructG <- function(knots,
                       scaleX,
                       pord) {
  if (pord == 2) {
    xi <- GrevillePoints(knots)
    G <- cbind(1, xi)
    if (scaleX) {
      q <- length(xi)
      xi_mn <- mean(xi)
      alpha <- normVec(xi-xi_mn)
      K <- matrix(data=c(1/sqrt(q),0,-xi_mn/alpha,1/alpha),nrow=2,ncol=2)
      G <- G %*% K
    }
  } else { # pord==1
    d <- attr(knots,"degree")
    q <- length(knots) - (d + 1)
    G <- matrix(data=1, nrow=q, ncol=1)
    if (scaleX) {
      G <- (1/sqrt(q))*G
    }
  }
  return(G)
}

#' Construct constraint matrix
#'
#' @param knots knot positions of B-spline basis.
#' @param pord order of the penalty matrix (pord=1 or 2).
#'
#' @return a q x q matrix of type spam
#'
#' @noRd
#' @keywords internal
constructCCt <- function(knots,
                         pord) {
  xmin <- attr(knots, which = "xmin")
  xmax <- attr(knots, which = "xmax")
  # if pord == 1 take point halfway, otherwise the
  # begin and endpoint of B-spline basis:
  if (pord == 1) {
    x <- 0.5 * (xmin + xmax)
  } else {
    x <- c(xmin, xmax)
  }
  Bx <- Bsplines(knots, x)
  CCt <- spam::crossprod.spam(Bx)
  return(CCt)
}

#' Construct list of Ginv matrices of splines
#'
#' @param q vector, with dimension of the B-spline basis used.
#' @param knots list with knot positions for each dimension
#' @param pord order of the penalty matrix (1 or 2).
#'
#' @return a list of symmetric matrices of length of vector q.
#'
#' @noRd
#' @keywords internal
constructGinvSplines <- function(q,
                                 knots,
                                 pord,
                                 scaleFactor) {
  ## dimension
  d <- length(q)
  lCCt <- lapply(X = seq_len(d),
                 FUN = function(i) { constructCCt(knots[[i]], pord)})
  CCt <- Reduce("kronecker", lCCt)
  lGinv <- list()
  for (i in seq_len(d)) {
    L <- list()
    for (j in seq_len(d)) {
      if (i == j) {
        L[[j]] <- scaleFactor[j]*constructPenalty(q[j], pord = pord)
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
#' @noRd
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
#' @noRd
#' @keywords internal
calcNomEffDim <- function(X,
                          Z,
                          dim.r,
                          term.labels.r) {
  if (is.null(Z)) return(NULL)
  p <- ncol(X)
  EDnom <- vector(length = length(dim.r))
  e <- cumsum(dim.r)
  s <- e - dim.r + 1
  # for each variance component in Z:
  for (i in 1:length(dim.r)) {
    ndx <- s[i]:e[i]
    Zi <- Z[, ndx, drop = FALSE]
    # if number of columns is high, use upper bound:
    if (dim.r[i] > 100 | nrow(X) > 10000) {
      rowSum <- spam::rowSums.spam(Zi)
      colSum <- spam::colSums.spam(Zi)
      nNonZeroCols <- sum(colSum != 0)
      if (all(rowSum == rowSum[1])) {
        EDnom[i] <- nNonZeroCols - 1
      } else {
        EDnom[i] <- nNonZeroCols
      }
    } else {
      XZ <- cbind(X, Zi)
      r <- qr(as.matrix(XZ))$rank
      if (r == p) {
        stop("Singularity problem with term ", term.labels.r[i],
             " in the random part of the model.\n",
             call. = FALSE)
      }
      EDnom[i] <- r - p
    }
  }
  return(EDnom)
}

#' Check variables in formula
#'
#' Check that all variables in a formula are present in the data. The variables
#' are converted to a factor in the data.frame as well.
#'
#' @importFrom utils hasName
#'
#' @noRd
#' @keywords internal
checkFormVars <- function(formula,
                          data,
                          naAllowed = TRUE) {
  if (!is.null(formula)) {
    ## Get variables in formula.
    formVars <- all.vars(formula)
    ## Get the name of formula as passed to the function.
    formName <- deparse(substitute(formula))
    missVars <- formVars[!hasName(data, formVars)]
    if (length(missVars) > 0) {
      stop("The following variables in the ", formName, " part of the model ",
           "are not in the data:\n", paste0(missVars, collapse = ", "), "\n",
           call. = FALSE)
    }
    if (length(formVars) > 0) {
      ## Check that variables have at least some non NA values.
      naVars <- formVars[sapply(X = formVars, FUN = function(formVar) {
        all(is.na(data[[formVar]]))
      })]
      if (length(naVars) > 0) {
        ## Get the name of formula as passed to the function.
        stop("The following variables in the ", formName, " part of the model ",
             "only have missing values:\n",
             paste0(naVars, collapse = ", "), "\n", call. = FALSE)
      }
      ## Check that variables have no NA values.
      if (!naAllowed) {
        ## NA allowed in response variable.
        respVar <- formVars[attr(terms(formula), "response")]
        nonRespVars <- setdiff(formVars, respVar)
        if (length(nonRespVars) > 0) {
          anyNaVars <- nonRespVars[sapply(X = nonRespVars,
                                          FUN = function(nonRespVar) {
                                            any(is.na(data[[nonRespVar]]))
                                          })]
          if (length(anyNaVars) > 0) {
            ## Get the name of formula as passed to the function.
            stop("The following variables in the ", formName,
                 " part of the model have missing values:\n",
                 paste0(anyNaVars, collapse = ", "), "\n", call. = FALSE)
          }
        }
      }
      for (formVar in formVars) {
        if (is.character(data[[formVar]])) {
          data[[formVar]] <- factor(data[[formVar]])
        }
      }
    }
  }
  return(data)
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
#' @noRd
#' @keywords internal
checkGroup <- function(random,
                       group,
                       data) {
  ## Only check if at least one of random and group is not NULL.
  if (!is.null(random) || !is.null(group)) {
    ## Check that all columns in group are in data.
    groupCols <- unlist(group)
    if (any(groupCols < 1) || any(groupCols > ncol(data))) {
      stop("All columns specified in group should be columns in data.\n",
           call. = FALSE)
    }
    if (is.null(random)) {
      grpVars <- NULL
    } else {
      ## Extract variables specified in grp() functions in random part.
      randTerms <- terms(random, specials = "grp")
      ## The next line seems redundant, but for some reason it reorders
      ## the terms/attributes in such a way that the order matches.
      ## The following code is based on a matching order.
      randTerms <- randTerms[seq_along(labels(randTerms))]
      grpPos <- attr(randTerms, which = "specials")$grp
      grpTerms <- sapply(X = grpPos, FUN = function(pos) {randTerms[pos]})
      grpVars <- sapply(X = grpTerms, FUN = all.vars)
      ## Check for variables in grp missing in group.
      if (length(grpVars) > 0) {
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
#' @noRd
#' @keywords internal
nameCoefs <- function(coefs,
                      desMat = NULL,
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
        attr(coefI, "termType") <- "(Intercept)"
      } else if (startsWith(x = labI, prefix = "lin(") ||
                 startsWith(x = labI, prefix = "s(") ||
                 startsWith(x = labI, prefix = "cf(")) {
        ## Spline and (cf) terms are just named 1...n.
        names(coefI) <- paste0(labI, "_", seq_along(coefI))
        attr(coefI, "termType") <- "fun" ## function
      } else if (!is.null(group) && labI %in% names(group)) {
        ## For group combine group name and column name.
        names(coefI) <- paste(labI, colnames(data)[group[[labI]]], sep = "_")
        attr(coefI, "termType") <- "grp"
      } else if (is.factor(data[[labI]])) {
        ## For fixed terms an extra 0 for the reference level has to be added.
        if (type == "fixed") {
          coefI <- c(0, coefI)
        }
        coefINames <- paste(labI, levels(data[[labI]]), sep = "_")
        ## Restrict to names in design matrix.
        if (!is.null(desMat)) {
          coefINames <-
            c(coefINames[1],
              coefINames[paste0(labI,
                                levels(data[[labI]])) %in% colnames(desMat)])
        }
        names(coefI) <- coefINames
        attr(coefI, "termType") <- "factor"
      } else {
        ## Numerical variable. Name equal to label.
        names(coefI) <- labI
        attr(coefI, "termType") <- "variable"
      }
    } else {
      ## Interaction term.
      ## Factors and non-factors in term have to be treated differently.
      isFactLab <- sapply(X = data[, labISplit], FUN = is.factor)
      labDatFact <- unique(data[names(isFactLab)[isFactLab]])
      labDatNonFact <- names(isFactLab)[!isFactLab]
      if (length(labDatFact) > 0 && length(labDatNonFact) > 0) {
        labDatNonFactDf <- as.data.frame(matrix(data = labDatNonFact,
                                                nrow = nrow(labDatFact),
                                                ncol = length(labDatNonFact),
                                                byrow = TRUE,
                                                dimnames = list(NULL, labDatNonFact)))
        labDat <- cbind(labDatFact, labDatNonFactDf)
        colnames(labDat) <- c(colnames(labDatFact), labDatNonFact)
        labDat <- labDat[labISplit]
      } else {
        labDat <- if (length(labDatFact) > 0) labDatFact else labDatNonFact
      }
      labDat <- lapply(X = seq_along(labDat), FUN = function(i) {
        if (isFactLab[i]) {
          paste0(labISplit[i], "_", labDat[[i]])
        } else {
          labDat[[i]]
        }
      })
      if (type == "fixed") {
        labDatAlt <- lapply(X = labDat, FUN = function(i) {
          gsub("_", "", i)
        })
        interAcLevsIAlt <- levels(interaction(labDatAlt, sep = ":",
                                              drop = TRUE))
        interAcLevsI <- levels(interaction(labDat, sep = ":", drop = TRUE))
        ## For fixed terms an extra 0 for the reference levels have to be added.
        coefI <- c(rep(0, times = length(interAcLevsIAlt) - length(coefI)),
                   coefI)
        ## Zeros correspond to missing columns from design matrix.
        ## Non-zeros to columns that are present in design matrix.
        names(coefI) <- c(interAcLevsI[!interAcLevsIAlt %in% colnames(desMat)],
                          interAcLevsI[interAcLevsIAlt %in% colnames(desMat)])
        ## Reorder.
        coefI <- coefI[interAcLevsI]
        names(coefI) <- levels(interaction(labDat, sep = ":", drop = TRUE))

      } else {
        names(coefI) <- levels(interaction(labDat, sep = ":", drop = FALSE))
      }
      attr(coefI, "termType") <- "interaction"
    }
    coefRes[[i]] <- coefI
  }
  return(coefRes)
}

#' Helper function to check Conditional Factor
#' Gives an error if syntax not correct, otherwise returns a boolean:
#' TRUE if a conditional factor is used, otherwise FALSE.
#'
#' @noRd
#' @keywords internal
checkConditionalFactor <- function(var,
                                   cond,
                                   level) {
  ## if conditional factor is defined
  vName <- deparse(substitute(var, env = parent.frame()))
  cName <- deparse(substitute(cond, env = parent.frame()))
  conditional <- FALSE
  if (vName != "NULL") {
    if (inherits(try(var, silent = TRUE), "try-error")) {
      stop(vName, " should be a variable in data.\n", call. = FALSE)
    }
  }
  if (cName != "NULL") {
    if (inherits(try(cond, silent = TRUE), "try-error") || is.null(cond)) {
      stop(cName, " should be a variable in data.\n", call. = FALSE)
    }
    if (!is.factor(cond)) {
      stop("cond should be a factor.\n", call. = FALSE)
    }
    if (is.null(level)) {
      stop("level should be defined.\n", call. = FALSE)
    }
    if (!(level %in% levels(cond))) {
      stop(level, " is not a level of ", cName, ".\n", call. = FALSE)
    }
    conditional <- TRUE
  }
  return(conditional)
}

#' Extend the Spam matrix with extra zero rows.
#'
#' @noRd
#' @keywords internal
extSpamMatrix <- function(X,
                          s,
                          N) {
  rPtr <- rep(0, N)
  rPtr[s] <- diff(X@rowpointers, differences = 1)
  rPtr <- c(0, rPtr)
  rPtr <- cumsum(rPtr) + 1
  ## change the slots of X.
  X@dimension[1] <- N
  X@rowpointers <- rPtr
  return(X)
}

#' Helper function for checking limits
#' The same check is done in the different spline functions so a helper
#' function is created to avoid duplicate code.
#'
#' @noRd
#' @keywords internal
checkLim <- function(lim,
                     limName,
                     x,
                     xName) {
  if (!is.numeric(lim)) {
    stop(limName, " should be numeric.\n", call. = FALSE)
  }
  if (length(lim) != 2) {
    stop(limName, " should be a vector of length two.\n", call. = FALSE)
  }
  if (lim[1] > range(x)[1] || lim[2] < range(x)[2]) {
    stop("All values of ", xName , " should be between the lower ",
         "and upper value of ", limName, ".\n", call. = FALSE)
  }
}

#' @noRd
#' @keywords internal
makeDesignTerm <- function(obj, newdat, term) {
  coefI <- obj$ndxCoefficients[[term]]
  type <- attr(coefI, which="termType")
  #cat("type: ", type, "\n")
  if (type == "factor") {
    lab <- names(coefI)
    df1 <- data.frame(lab, ndx = as.numeric(coefI))
    x <- paste0(term, "_", newdat[[term]])
    chk <- x %in% lab
    if (sum(!chk) > 0) {
      lev <- paste(unique(x[!chk]), collapse=',')
      err <- paste("levels", lev, "not defined")
      stop(err)
    }
    df2 <- data.frame(lab=x, nr=c(1:length(x)))
    df <- merge(df1, df2)
    df <- df[df$ndx!=0, ]
    l <- list(i=df$nr,j=df$ndx,v = rep(1, nrow(df)))
    M <- spam::spam(l, nrow=length(x), ncol=nrow(obj$C))
  } else if (type == "variable") {
    v <- newdat[[term]]
    ndx <- as.numeric(coefI)
    l <- list(i=1:length(v), j=rep(ndx,length(v)), v=v)
    M <- spam::spam(l, nrow=length(v), ncol=nrow(obj$C))

  } else {
    str <- paste("Make predictions for type", type, "not implemented yet\n" )
    stop(str)
  }
  return(M)
}


