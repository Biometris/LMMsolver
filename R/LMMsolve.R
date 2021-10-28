#' Solve Linear Mixed Models
#'
#' Solve Linear Mixed Models.
#'
#' @param fixed A formula for the fixed part of the model. Should be of the
#' form "response ~ pred"
#' @param random A formula for the random part of the model. Should be of the
#' form "~ pred".
#' @param spline A formula for the spline part of the model. Should be of the
#' form "~ spl1D()", ~ spl2D()" or "~spl3D()".
#' @param group A named list where each component is a numeric vector
#' specifying contiguous fields in data that are to be considered as a
#' single term.
#' @param ginverse A named list with each component a symmetric matrix, the
#' inverse of the vcov matrix of a random term in the model.
#' @param data A data.frame containing the modeling data.
#' @param residual A formula for the residual part of the model. Should be of
#' the form "~ pred".
#' @param tolerance A numerical value. The convergence tolerance for the
#' modified Henderson algorithm to estimate the variance components.
#' @param trace Should the progress of the algorithm be printed? Default
#' \code{trace = FALSE}.
#' @param maxit A numerical value. The maximum number of iterations for the
#' algorithm. Default \code{maxit = 250}.
#' @param theta initial values for penalty or precision parameters. Default
#' \code{NULL}, all precision parameters set equal to 1.
#'
#' @return An object of class \code{LMMsolve} representing the fitted model.
#' See \code{\link{LMMsolveObject}} for a full description of the components in
#' this object.
#'
#' @examples
#' ## Fit models on john.alpha data from agridat package.
#' data(john.alpha, package = "agridat")
#'
#' ## Fit simple model with only fixed effects.
#' LMM1 <- LMMsolve(fixed = yield ~ rep + gen,
#'                 data = john.alpha)
#'
#' ## Fit the same model with genotype as random effect.
#' LMM1_rand <- LMMsolve(fixed = yield ~ rep,
#'                      random = ~gen,
#'                      data = john.alpha)
#'
#' ## Fit the model with a 1-dimensional spline at the plot level.
#' LMM1_spline <- LMMsolve(fixed = yield ~ rep + gen,
#'                        spline = ~spl1D(x = plot, nseg = 20),
#'                        data = john.alpha)
#'
#' ## Fit models on multipop data included in the package.
#' data(multipop)
#'
#' ## The residual variances for the two populations can be different.
#' ## Allow for heterogeneous residual variances using the residual argument.
#' LMM2 <- LMMsolve(fixed = pheno ~ cross,
#'                 residual = ~cross,
#'                 data = multipop)
#'
#' ## QTL-probabilities are defined by the columns pA, pB, pC.
#' ## They can be included in the random part of the model by specifying the
#' ## group argument and using grp() in the random part.
#'
#' # Define groups by specifying columns in data corresponding to groups in a list.
#' # Name used in grp() should match names specified in list.
#' lGrp <- list(QTL = 3:5)
#' LMM2_group <- LMMsolve(fixed = pheno ~ cross,
#'                       group = lGrp,
#'                       random = ~grp(QTL),
#'                       residual = ~cross,
#'                       data = multipop)
#'
#' @seealso \code{\link{LMMsolveObject}}, \code{\link{spl1D}},
#' \code{\link{spl2D}}, \code{\link{spl3D}}
#'
#' @importFrom stats model.frame terms model.matrix contrasts as.formula
#' terms.formula aggregate
#'
#' @export
LMMsolve <- function(fixed,
                     random = NULL,
                     spline = NULL,
                     group = NULL,
                     ginverse = NULL,
                     data,
                     residual = NULL,
                     tolerance = 1.0e-6,
                     trace = FALSE,
                     maxit = 250,
                     theta = NULL) {
  ## Input checks.
  if (!inherits(data, "data.frame")) {
    stop("data should be a data.frame.\n")
  }
  if (!inherits(fixed, "formula") || length(terms(fixed)) != 3) {
    stop("fixed should be a formula of the form \"resp ~ pred\".\n")
  }
  if (!is.null(random) &&
      (!inherits(random, "formula") || length(terms(random)) != 2)) {
    stop("random should be a formula of the form \" ~ pred\".\n")
  }
  if (!is.null(spline) &&
      (!inherits(spline, "formula") || length(terms(spline)) != 2 ||
       ## Spline formula should consist of splxD() and nothing else.
       length(unlist(attr(terms(spline, specials = c("spl1D", "spl2D", "spl3D")),
                          "specials"))) != 1)) {
    stop("spline should be a formula of form \"~ spl1D()\", \"~ spl2D()\" ",
         "or \"~spl3D()\".\n")
  }
  if (!is.null(ginverse) &&
      (!is.list(ginverse) ||
       length(names(ginverse)) == 0 ||
       (!all(sapply(X = ginverse, FUN = function(x) {
         is.matrix(x) && isSymmetric(x)}))))) {
    stop("ginverse should be a named list of symmetric matrices.\n")
  }
  if (!is.null(residual) &&
      (!inherits(residual, "formula") || length(terms(residual)) != 2)) {
    stop("residual should be a formula of the form \" ~ pred\".\n")
  }
  if (!is.numeric(tolerance) || length(tolerance) > 1 || tolerance < 0) {
    stop("tolerance should be a positive numerical value.")
  }
  if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 0) {
    stop("maxit should be a positive numerical value.")
  }
  ## Drop unused factor levels from data.
  data <- droplevels(data)
  ## Check that all variables used in formulas are in data.
  random <- checkGroup(random, group)
  checkFormVars(fixed, data)
  checkFormVars(random, data)
  checkFormVars(residual, data)
  ## Remove NA for response variable from data.
  respVar <- all.vars(fixed)[attr(terms(fixed), "response")]
  respVarNA <- is.na(data[[respVar]])
  if (sum(respVarNA) > 0) {
    warning(sum(respVarNA), " observations removed with missing value for ",
            respVar, ".\n", call. = FALSE)
    data <- data[!respVarNA, ]
  }
  ## Make random part.
  if (!is.null(random)) {
    mf <- model.frame(random, data, drop.unused.levels = TRUE, na.action = NULL)
    mt <- terms(mf)
    f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
    Z1 <- model.matrix(mt, data = mf,
                       contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                              FUN = contrasts,
                                              contrasts = FALSE))
    dim1.r <- table(attr(Z1, "assign"))[-1]
    term1.labels.r <- attr(mt, "term.labels")
    ## Number of variance parameters (see Gilmour 1995) for each variance component
    varPar1 <- rep(1, length(dim1.r))
    Z1 <- Z1[, -1]
  } else {
    dim1.r <- NULL
    term1.labels.r <- NULL
    Z1 <- NULL
    varPar1 <- NULL
  }
  if (!is.null(group)) {
    ndx <- unlist(group)
    dim2.r <- sapply(X = group, FUN = length)
    term2.labels.r <- names(group)
    varPar2 <- rep(1, length(dim2.r))
    Z2 <- as.matrix(data[, ndx])
  } else {
    dim2.r <- NULL
    term2.labels.r <- NULL
    Z2 <- NULL
    varPar2 <- NULL
  }
  if (!(is.null(random) & is.null(group))) {
    Z <- cbind(Z1, Z2)
    dim.r <- c(dim1.r, dim2.r)
    term.labels.r <- c(term1.labels.r, term2.labels.r)
    varPar <- c(varPar1, varPar2)
    e <- cumsum(dim.r)
    s <- e - dim.r + 1
    lGinv <- list()
    for (i in 1:length(dim.r)) {
      tmp <- rep(0, sum(dim.r))
      tmp[s[i]:e[i]] <- 1
      lGinv[[i]] <- spam::diag.spam(tmp)
    }
    names(lGinv) <- term.labels.r
  } else {
    Z <- NULL
    lGinv <- NULL
    dim.r <- NULL
    term.labels.r <- NULL
    varPar <- NULL
  }
  ## Replace identity matrix with ginverse for random terms.
  if (!is.null(ginverse)) {
    ## Convert to spam.
    e <- cumsum(dim.r)
    s <- e - dim.r + 1
    for (i in seq_along(ginverse)) {
      ginvVar <- names(ginverse)[i]
      ginvMat <- ginverse[[i]]
      k <- which(term.labels.r == ginvVar)
      if (length(k) == 0) {
        stop("ginverse element ", ginvVar, " not defined in random part.\n")
      }
      if (dim.r[k] != nrow(ginvMat)) {
        stop("Dimensions of ", ginvVar, " should match number of levels ",
             "for corresponding factor in data.\n")
      }
      if (!isTRUE(all(rownames(ginvMat) == levels(data[[ginvVar]]))) ||
          !isTRUE(all(colnames(ginvMat) == levels(data[[ginvVar]])))) {
        stop("Row and column names of ", ginvVar, " should match levels ",
             "of corresponding factor in data.\n")
      }
      ndx <- s[k]:e[k]
      ## as.spam drops row and column names, so only converting to spam
      ## after all checking is done.
      lGinv[[k]][ndx, ndx] <- spam::as.spam(ginvMat)
    }
  }
  ## Make fixed part.
  mf <- model.frame(fixed, data, drop.unused.levels = TRUE)
  mt <- terms(mf)
  f.terms <- all.vars(mt)[attr(mt, "dataClasses") == "factor"]
  X = model.matrix(mt, data = mf,
                   contrasts.arg = lapply(X = mf[, f.terms, drop = FALSE],
                                          FUN = contrasts, contrasts = TRUE))
  dim.f <- as.numeric(table(attr(X, "assign")))
  term.labels.f <- attr(mt, "term.labels")
  # calculate NomEff dimension for non-spline part
  NomEffDimRan <- calcNomEffDim(X, Z, dim.r)
  ## Add spline part.
  splRes <- NULL
  if (!is.null(spline)) {
    tf <- terms(spline, specials = c("spl1D", "spl2D", "spl3D"))
    terms <- attr(tf, "term.labels")
    splRes <- eval(parse(text = terms), envir = data, enclos = parent.frame())
    ## Add to design matrix fixed effect X.
    X <- cbind(X, splRes$X)
    ## Add to design matrix random effect Z.
    Z <- cbind(Z, splRes$Z)
    ## Expand matrices Ginv to the updated Z.
    lGinv <- expandGinv(lGinv, splRes$lGinv)
    ## A splxD model has x parameters.
    varPar <- c(varPar, length(splRes$lGinv))
    ## Add dims.
    dim.f <- c(dim.f, splRes$dim.f)
    dim.r <- c(dim.r, splRes$dim.r)
    ## Add nominal ED
    NomEffDimRan <- c(NomEffDimRan, splRes$EDnom)
    ## Add labels.
    term.labels.f <- c(term.labels.f, splRes$term.labels.f)
    term.labels.r <- c(term.labels.r, splRes$term.labels.r)
  }
  ## Add intercept.
  if (attr(mt, "intercept") == 1) {
    term.labels.f <- c("(Intercept)", term.labels.f)
  }
  ## Make residual part.
  if (!is.null(residual)) {
    residVar <- all.vars(residual)
    lRinv <- makeRlist(df = data, column = residVar)
  } else {
    lRinv <- list(residual = spam::diag.spam(1, nrow(data)))
    attr(lRinv, "cnt") <- nrow(data)
  }
  y <- mf[, 1]
  obj <- sparseMixedModels(y = y, X = X, Z = Z, lGinv = lGinv, lRinv = lRinv,
                           tolerance = tolerance, trace = trace, maxit = maxit,
                           theta = theta)
  ## Add names to coefficients.
  ## Fixed terms.
  ef <- cumsum(dim.f)
  sf <- ef - dim.f + 1
  coefF <- vector(mode = "list", length = length(dim.f))
  for (i in seq_along(coefF)) {
    coefFi <- obj$a[sf[i]:ef[i]]
    labFi <- term.labels.f[i]
    if (labFi == "(Intercept)") {
      names(coefFi) <- "(Intercept)"
      ## For fixed terms an extra 0 for the reference level has to be added.
    } else if (labFi == "splF") {
      ## Spline terms are just named 1...n.
      names(coefFi) <- paste0("splF_", seq_along(coefFi))
    } else {
      coefFi <- c(0, coefFi)
      names(coefFi) <- paste(labFi, levels(data[[labFi]]) , sep = "_")
    }
    coefF[[i]] <- coefFi
  }
  ## Similar for random terms.
  er <- sum(dim.f) + cumsum(dim.r)
  sr <- er - dim.r + 1
  coefR <- vector(mode = "list", length = length(dim.r))
  for (i in seq_along(coefR)) {
    coefRi <- obj$a[sr[i]:er[i]]
    labRi <- term.labels.r[i]
    if (labRi == "splR") {
      ## Spline terms are just named 1...n.
      names(coefRi) <- paste0("splR_", seq_along(coefRi))
    } else if (labRi %in% names(group)) {
      ## For group combine group name and column name.
      names(coefRi) <- paste(labRi, colnames(data)[group[[labRi]]], sep = "_")
    } else {
      names(coefRi) <- paste(labRi, levels(data[[labRi]]), sep = "_")
    }
    coefR[[i]] <- coefRi
  }
  ## Combine result for fixed and random terms.
  coefTot <- c(coefF, coefR)
  names(coefTot) <- c(term.labels.f, term.labels.r)
  ## Extract effective dimensions from fitted model.
  EffDimRes <- attributes(lRinv)$cnt
  EffDimNamesRes <- attributes(lRinv)$names
  NomEffDim <- c(NomEffDimRan, EffDimRes)
  # Calc upper bound for nominal effective dimension:
  N <- nrow(X)
  p <- ncol(X)
  NomEffDim <- pmin(NomEffDim, N - p)
  ## Make ED table for fixed effects.
  EDdf1 <- data.frame(Term = term.labels.f,
                      Effective = dim.f,
                      Model = dim.f,
                      Nominal = dim.f,
                      Ratio = rep(1, length(dim.f)),
                      Penalty = rep(0, length(dim.f)),
                      VarComp = rep(NA, length(dim.f)))
  ## Make ED table for random effects.
  EDdf2 <- data.frame(Term = obj$EDnames,
                      Effective = obj$ED,
                      Model = c(rep(dim.r, varPar), EffDimRes),
                      Nominal = NomEffDim,
                      Ratio = obj$ED / NomEffDim,
                      Penalty = obj$theta,
                      VarComp = c(rep(term.labels.r, varPar), EffDimNamesRes))
  ## Make full ED table.
  EDdf <- rbind(EDdf1, EDdf2)
  EDdf <- EDdf[-which(colnames(EDdf) == "VarComp")]
  rownames(EDdf) <- NULL
  ## Make variance table.
  varComp <- factor(EDdf2[["VarComp"]], levels = unique(EDdf2[["VarComp"]]))
  VarDf <- aggregate(x = EDdf2$Penalty, by = list(VarComp = varComp), FUN = sum)
  VarDf[["Variance"]] = 1 / VarDf[["x"]]
  VarDf <- VarDf[c("VarComp", "Variance")]
  ## Compute REML constant.
  constantREML <- -0.5 * log(2 * pi) * (length(y) - sum(dim.f))
  dim <- c(dim.f, dim.r)
  return(LMMsolveObject(logL = obj$logL,
                        sigma2e = obj$sigma2e,
                        tau2e = obj$tau2e,
                        EDdf = EDdf,
                        varPar = varPar,
                        VarDf = VarDf,
                        theta = obj$theta,
                        coefficients = coefTot,
                        yhat = obj$yhat,
                        residuals = obj$residuals,
                        nIter = obj$nIter,
                        C = obj$C,
                        constantREML = constantREML,
                        dim = dim,
                        Nres = length(lRinv),
                        term.labels.f = term.labels.f,
                        term.labels.r = term.labels.r,
                        splRes = splRes))
}
