
#' define cf() function, to allow for variances conditional on factor:
#' @noRd
cf <- function(var, cond, level) {
  vName <- deparse(substitute(var))
  cName <- deparse(substitute(cond))

  f <- as.formula(paste0("~-1+",vName))

  mf <- model.frame(f, dat, drop.unused.levels = TRUE, na.action = NULL)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
  Z <- Matrix::sparse.model.matrix(mt, data = mf,
                                    contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                                           FUN = contrasts,
                                                           contrasts = FALSE))
  Z <- Z * (cond == level)
  ndx <- which(spam::colSums(Z)!=0)
  Z <- Z[, ndx]
  Z <- spam::as.spam.dgCMatrix(Z)
  termlabel <- paste0("cf(",cName,", ",level,")_",vName)
  cflabel <- paste(termlabel, levels(var)[ndx],sep="_")
  L <- list(Z=Z, termlabel = termlabel, cflabel=cflabel, dim.r=ncol(Z))
  return(L)
}

#' Function to analyse cf terms in the random part.
#' @keywords internal
condFactor <- function(random, data) {
  tf <- terms.formula(random, specials = c("cf"))
  f <- attr(tf, "term.labels")
  ndxAt <- attr(tf,"specials")$cf
  if (is.null(ndxAt)) return(NULL)

  g <- f[ndxAt]
  Nterms <- length(g)
  dim.r <- NULL
  Z <- NULL
  termlabel <- NULL
  cflabel <- list()
  for (i in 1:Nterms) {
    L <- eval(parse(text = g[i]), envir = data, enclos = parent.frame())
    Z <- cbind.spam(Z, L$Z)
    termlabel <- c(termlabel, L$termlabel)
    cflabel[[L$termlabel]] <- L$cflabel
    dim.r <- c(dim.r, L$dim.r)
  }
  random <- tf[-ndxAt]
  return(list(Z=Z, term.labels.r=termlabel, coefLabels=cflabel, dim.r=dim.r, random=random))
}


