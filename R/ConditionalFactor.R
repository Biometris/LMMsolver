
#' define cf() function, to allow for variances conditional on factor:
#'
#' @noRd
cf <- function(var,
               cond,
               level) {
  vName <- deparse(substitute(var))
  cName <- deparse(substitute(cond))
  ## check (syntax) conditional factor.
  checkConditionalFactor(var, cond, level)
  ndx <- cond == level
  var <- droplevels(var[ndx])
  df <- data.frame(x = var)
  Z <- Matrix::sparse.model.matrix(~-1 + x, data = df)
  Z <- spam::as.spam.dgCMatrix(Z)
  s <- which(ndx)
  Z <- extSpamMatrix(X = Z, s = s, N = length(ndx))
  termlabel <- paste0("cf(", cName, ", ", level, ")_", vName)
  cflabel <- paste(termlabel, levels(var), sep = "_")
  L <- list(Z = Z, termlabel = termlabel, cflabel = cflabel, dim.r = ncol(Z))
  return(L)
}

#' Function to analyse cf terms in the random part.
#'
#' @noRd
#' @keywords internal
condFactor <- function(random,
                       data) {
  if (is.null(random)) return(NULL)
  tf <- terms.formula(random, specials = "cf")
  f <- attr(tf, "term.labels")
  ndxAt <- attr(tf, "specials")$cf
  if (is.null(ndxAt)) return(NULL)
  g <- f[ndxAt]
  Nterms <- length(g)
  dim.r <- NULL
  Z <- NULL
  termlabel <- NULL
  cflabel <- list()
  for (gTerm in g) {
    L <- eval(parse(text = gTerm), envir = data, enclos = parent.frame())
    Z <- spam::cbind.spam(Z, L$Z)
    termlabel <- c(termlabel, L$termlabel)
    cflabel[[L$termlabel]] <- L$cflabel
    dim.r <- c(dim.r, L$dim.r)
  }
  random <- tf[-ndxAt]
  return(list(Z = Z, term.labels.r = termlabel, coefLabels = cflabel,
              dim.r = dim.r, random = random))
}


