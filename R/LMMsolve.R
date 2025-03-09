#' Solve Linear Mixed Models
#'
#' Solve Linear Mixed Models using REML.
#'
#' A Linear Mixed Model (LMM) has the form
#' \deqn{y = X \beta + Z u + e, u \sim N(0,G), e \sim N(0,R)} where
#' \eqn{y} is a vector of observations, \eqn{\beta} is a vector with the fixed
#' effects, \eqn{u} is a vector with the random effects, and \eqn{e} a vector of
#' random residuals. \eqn{X} and \eqn{Z} are design matrices.
#'
#' LMMsolve can fit models where the matrices \eqn{G^{-1}} and \eqn{R^{-1}} are
#' a linear combination of precision matrices \eqn{Q_{G,i}} and \eqn{Q_{R,i}}:
#' \deqn{G^{-1} = \sum_{i} \psi_i Q_{G,i} \;, R^{-1} = \sum_{i} \phi_i Q_{R,i}}
#' where the precision parameters \eqn{\psi_i} and \eqn{\phi_i} are estimated
#' using REML. For most standard mixed models \eqn{1/{\psi_i}} are the variance
#' components and \eqn{1/{\phi_i}} the residual variances. We use a formulation
#' in terms of precision parameters to allow for non-standard mixed models using
#' tensor product splines.
#'
#' @param fixed A formula for the fixed part of the model. Should be of the
#' form "response ~ pred"
#' @param random A formula for the random part of the model. Should be of the
#' form "~ pred".
#' @param spline A formula for the spline part of the model. Should be of the
#' form "~ spl1D()", ~ spl2D()" or "~spl3D()". Generalized Additive Models (GAMs) can
#' also be used, for example "~ spl1D() + spl2D()"
#' @param group A named list where each component is a numeric vector
#' specifying contiguous fields in data that are to be considered as a
#' single term.
#' @param ginverse A named list with each component a symmetric matrix, the
#' precision matrix of a corresponding random term in the model. The row and
#' column order of the precision matrices should match the order of the
#' levels of the corresponding factor in the data.
#' @param weights A character string identifying the column
#' of data to use as relative weights in the fit. Default value NULL, weights are
#' all equal to one.
#' @param data A data.frame containing the modeling data.
#' @param residual A formula for the residual part of the model. Should be of
#' the form "~ pred".
#' @param family An object of class \code{family} or \code{familyLMMsolver}
#' specifying the distribution and link function. See class \code{\link{family}} and
#' and \code{\link{multinomial}} for details.
#' @param offset An a priori known component to be included in the linear
#' predictor during fitting. \code{Offset} be a numeric vector, or a character
#' string identifying the column of data. Default \code{offset = 0}.
#' @param tolerance A numerical value. The convergence tolerance for the
#' modified Henderson algorithm to estimate the variance components.
#' @param trace Should the progress of the algorithm be printed? Default
#' \code{trace = FALSE}.
#' @param maxit A numerical value. The maximum number of iterations for the
#' algorithm. Default \code{maxit = 250}.
#' @param theta initial values for penalty or precision parameters. Default
#' \code{NULL}, all precision parameters set equal to 1.
#' @param grpTheta a vector to give components the same penalty. Default
#' \code{NULL}, all components have a separate penalty.
#'
#' @returns An object of class \code{LMMsolve} representing the fitted model.
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
#' terms.formula aggregate model.response var formula gaussian
#'
#' @export
LMMsolve <- function(fixed,
                     random = NULL,
                     spline = NULL,
                     group = NULL,
                     ginverse = NULL,
                     weights = NULL,
                     data,
                     residual = NULL,
                     family = gaussian(),
                     offset = 0,
                     tolerance = 1.0e-6,
                     trace = FALSE,
                     maxit = 250,
                     theta = NULL,
                     grpTheta = NULL) {
  ## Input checks, and return numeric value for offset
  offset <- chkInputLMMsolve(fixed = fixed, random = random,
              data = data, ginverse = ginverse,
              residual = residual, tolerance = tolerance,
              maxit = maxit, grpTheta = grpTheta, offset = offset,
              family = family)

  ## Check that all variables used in fixed formula are in data.
  data <- checkFormVars(fixed, data, naAllowed = FALSE)

  mult_col_response <- FALSE
  if (family$family %in% c("multinomial", "binomial", "quasibinomial")) {
    mf <- model.frame(fixed, data, drop.unused.levels = TRUE, na.action = NULL)
    respNames <- colnames(mf[[1]])
    chkFac <- sapply(respNames, function(x) {is.factor(data[[x]])})
    if (any(chkFac)) {
      str <- paste("response", names(mf)[1], "should be numeric.")
      stop(str)
    }
    YY <- model.response(mf, type = "any")
    if (family$family == "multinomial") {
      if (!is.null(weights)) {
        str <- paste("family", family$family,
                     ": weights option cannot be used.")
        stop(str)
      }
      checkMultiResponse(YY, family)
      respVar <- colnames(YY)
      respVarNA <- rep(FALSE, nrow(data))
      mult_col_response <- TRUE
      nCat <- ncol(YY) - 1
    } else { # binomial, quasibinomial
      if (inherits(YY, "matrix")) {
        mult_col_response <- TRUE
        checkMultiResponse(YY, family)
        if (!is.null(weights)) {
          str <- paste("family", family$family,
                       ": weights cannot be used in cbind format.")
          stop(str)
        }
        cntYY <- rowSums(YY)
        y_proportion <- YY[,1]/cntYY
        weights <- "weights"
        data$weights <- cntYY
        data$y_proportion <- y_proportion
        respVar <- "y_proportion"
        respVarNA <- rep(FALSE, nrow(data))
        fixed_txt <- as.character(fixed)
        fixed <- as.formula(paste0("y_proportion ~ ", fixed_txt[3]))
      } else {
        if (inherits(YY, "factor")) {
          respVar <- all.vars(fixed)[attr(terms(fixed), "response")]
          Factor <- data[[respVar]]
          levels <- levels(Factor)
          sc <- 1.0*(Factor == levels[2])
          data[[respVar]] <- sc
        }
      }
    }
  }
  if (!mult_col_response) {
    respVar <- all.vars(fixed)[attr(terms(fixed), "response")]
    ## Remove NA for response variable from data.
    respVarNA <- is.na(data[[respVar]])
    if (sum(respVarNA) > 0) {
      warning(sum(respVarNA), " observations removed with missing value for ",
              respVar, ".\n", call. = FALSE)
      data <- data[!respVarNA, ]
    }
  }
  ## Check for weights should be done after removing observations with NA for
  ## response var. This prevents error messages when both response var and
  ## weight are NA.
  w <- getWeights(weights, data)
  ## Remove observations with zero weights
  weightsZero <- w == 0
  if (sum(weightsZero) > 0) {
    # warning(sum(weightsZero), " observations removed with zero weights \n ",
    #        call. = FALSE)
    data <- data[!weightsZero, ]
    ## remove missing values for weight (default w=1).
    w <- w[!weightsZero]
  }
  ## Drop unused factor levels from data.
  data <- droplevels(data)

  ## Check random term for conditional factor
  condFactor <- condFactor(random, data)
  if (!is.null(condFactor)) {
    random <- condFactor$random
  }
  ## Check that all variables used in formulas are in data.
  chkGroup <- checkGroup(random, group, data)
  random <- chkGroup$random
  group <- chkGroup$group
  data <- checkFormVars(random, data)
  data <- checkFormVars(residual, data)

  ## construct random part of model (excluding splines)
  ranPart <- constructRandom(random, group, condFactor, data)
  Z <- ranPart$Z
  lGinv <- ranPart$lGinv
  dim.r <- ranPart$dim.r
  term.labels.r <- ranPart$term.labels.r
  scFactor <- ranPart$scFactor
  varPar <- ranPart$varPar
  nNonSplinesRandom <- ranPart$nNonSplinesRandom

  lGinv <- constructGinv(ginverse, lGinv, dim.r, term.labels.r)

  ## construct Fixed part (excluding splines)
  X <- constructFixed(fixed, data)
  dim.f <- attr(X, which="dim.f")
  term.labels.f <- attr(X, which="term.labels.f")

  ## Add spline part.
  splResList <- NULL
  NomEffDimRan <- NULL
  if (!is.null(spline)) {
    # check correct formula splines
    chkSplinesFormula(spline)
    splTrms <- terms(spline, specials = c("spl1D", "spl2D", "spl3D"))
    splSpec <- attr(splTrms, "specials")
    splTerms <- labels(splTrms)
    Nterms <- length(splTerms)
    splResList <- list()
    for (i in seq_len(Nterms)) {
      splRes <- eval(parse(text = splTerms[i]), envir = data, enclos = parent.frame())
      ## Multiple 1D gam models should have unique x variables.
      if (!is.null(term.labels.f) && !is.null(splRes$term.labels.f) &&
          splRes$term.labels.f %in% term.labels.f) {
        stop("x variables in 1D splines should be unique.\n")
      }
      splResList[[i]] <- splRes
      ## Add to design matrix fixed effect X.
      X <- cbind(X, splRes$X)
      ## Add to design matrix random effect Z.
      Z <- spam::cbind.spam(Z, splRes$Z)
      ## Expand matrices Ginv to the updated Z.
      lGinv <- expandGinv(lGinv, splRes$lGinv)
      ## A splxD model has x parameters.
      varPar <- c(varPar, length(splRes$lGinv))
      ## Add dims.
      dim.f <- c(dim.f, splRes$dim.f)
      dim.r <- c(dim.r, splRes$dim.r)
      ## Add nominal ED.
      if (family$family == "multinomial") {
        NomEffDimRan <- c(NomEffDimRan, splRes$EDnom*nCat)
      } else {
        NomEffDimRan <- c(NomEffDimRan, splRes$EDnom)
      }
      ## Add labels.
      term.labels.f <- c(term.labels.f, splRes$term.labels.f)
      term.labels.r <- c(term.labels.r, splRes$term.labels.r)
      scFactor <- c(scFactor, splRes$scaleFactor)
    }
    # if splines of different dimensions are combined:
    if (length(splSpec[!sapply(X = splSpec, FUN = is.null)]) > 1) {
       #check whether variables are uniquely defined
        nameVar <- unlist(lapply(splResList, function(z) {names(z$x)}))
        if (length(nameVar)!=length(unique(nameVar))) {
        stop("variables in splines should be unique.\n")
      }
    }
  }

  ## get the response variable, after all filtering has been done
  if (family$family != "multinomial") {
    y <- data[[respVar]]  # vector
    ## check whether the variance for response is not zero.
    chkResponse(y, residual, data)
  } else {
    y <- YY # matrix
  }

  ## Fit the model
  LMMobj <- fitLMM(y = y, X = X, Z = Z, w = w, lGinv = lGinv,
              tolerance = tolerance, trace = trace, maxit = maxit,
              theta = theta, grpTheta = grpTheta,
              family = family, offset = offset,
              dim.f = dim.f, dim.r = dim.r,
              term.labels.f = term.labels.f, term.labels.r = term.labels.r,
              respVar = respVar,
              NomEffDimRan = NomEffDimRan, varPar = varPar,
              splResList = splResList, residual = residual,
              group = group, nNonSplinesRandom = nNonSplinesRandom,
              scFactor = scFactor, data = data)
  return(LMMobj)
}
