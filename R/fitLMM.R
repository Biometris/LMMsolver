fitLMM <- function(y, X, Z, w, lGinv, tolerance, trace, maxit,
              theta, grpTheta, family, offset, dim.f, dim.r,
              term.labels.f, term.labels.r, respVar, NomEffDimRan,
              varPar, splResList, residual, group,
              nNonSplinesRandom, scFactor, data) {
  ## Convert to spam matrix and cleanup
  Xs <- spam::as.spam.dgCMatrix(X)
  Xs <- spam::cleanup(Xs)

  ## construct inverse of residual matrix R.
  lRinv <- constructRinv(df = data, residual = residual, weights = w)
  nRes <- length(lRinv)
  scFactor <- c(scFactor, rep(1, nRes))
  ## set theta
  if (!is.null(theta)) {
    if (length(theta) != length(scFactor)) {
      stop("Argument theta has wrong length \n")
    }
    theta <- theta / scFactor
  } else {
    theta <- 1 / scFactor
  }
  ## set grpTheta
  if (!is.null(grpTheta)) {
    if (length(grpTheta) != length(scFactor)) {
      stop("Argument grpTheta has wrong length \n")
    }
  } else {
    grpTheta <- c(1:length(scFactor))
  }
  if (family$family == "gaussian") {
    obj <- sparseMixedModels(y = y, X = Xs, Z = Z, lGinv = lGinv, lRinv = lRinv,
                             tolerance = tolerance, trace = trace, maxit = maxit,
                             theta = theta, grpTheta = grpTheta)
    dev.residuals <- family$dev.resids(y, obj$yhat, w)
    deviance <- sum(dev.residuals)
  } else if (family$family == "multinomial") {
    YY <- y
    n <- nrow(YY)
    nCat <- length(respVar) - 1
    sY <- rowSums(YY)
    YY <- YY/sY
    YY <- YY[,-(nCat+1)] # remove last column
    y <- as.vector(t(YY))
    Xs <- Xs %x% spam::diag.spam(nCat)
    Z <- Z %x% spam::diag.spam(nCat)
    lGinv <- lapply(lGinv, FUN = function(x) {x %x% spam::diag.spam(nCat)})

    # remark: We can add some prior information to initial
    # values based on the scores:
    eta <- rep(0, n*nCat)

    # Initialize the sparse block structure for W and Dinv
    W <- diag.spam(1, n) %x% spam(x=1, nrow=nCat, ncol=nCat)
    Dinv <- diag.spam(1, n) %x% spam(x=1, nrow=nCat, ncol=nCat)

    theta <- c(1,1)                # check this!
    fixedTheta <- c(FALSE, TRUE)   # check this!
    trace_GLMM <- NULL
    for (i in 1:maxit) {
      Eta <- matrix(data=eta, ncol=nCat, nrow=n, byrow=TRUE)
      Pi <- t(apply(X=Eta, MARGIN=1, FUN=h))
      pi <- as.vector(t(Pi))

      D_list <- lapply(1:n, FUN =
                         function(r) {
                           return(Jacobian(Eta[r,]))
                         } )
      Dinv_list <- lapply(1:n, FUN= function(r) { solve(D_list[[r]])})

      W_list <- lapply(1:n, FUN =
                         function(r) {
                           Sigma <- calcSigma(Pi[r,])
                           D_list[[r]] %*% solve(Sigma) %*% t(D_list[[r]])
                         } )

      # Store in the defined sparse format for W and Dinv:
      Dinv@entries <- unlist(Dinv_list)
      W@entries <- unlist(W_list)

      z <- eta + Dinv %*% (y - pi)
      lRinv <- list(W)
      attr(lRinv, "cnt") <- n # correct?
      obj <- sparseMixedModels(z, X = Xs, Z = Z,
                               lGinv = lGinv, lRinv = lRinv, trace=TRUE,
                               fixedTheta = fixedTheta, theta = theta)
      eta_old  <- eta
      eta <- obj$yhat
      theta <- obj$theta
      tol <- sum((eta - eta_old)^2)/sum(eta^2)

      aux.df <- data.frame(itOuter = rep(i, obj$nIter),
                           tol=c(rep(NA,obj$nIter-1), tol))
      trace_GLMM <- rbind(trace_GLMM, cbind(aux.df, obj$trace))

      cat("iter ", i, "     ", tol, "    ", obj$logL, "   ", obj$ED, "\n")
      if (tol < 1.0e-10) {
        cat("convergence after", i, "iterations\n")
        break;
      }
    }
    dev.residuals <- NULL # check how to calculate this!
    deviance <- NULL      # check how to calculate this!

    #stop("Still Working on multinomial() implementation!")
  } else {
    ## MB, 23 jan 2023
    ## binomial needs global weights
    weights <- w
    nobs <- length(y)
    mustart <- etastart <- NULL
    eval(family$initialize)
    mu <- mustart
    eta <- family$linkfun(mustart)
    nNonRes <- length(theta) - nRes
    fixedTheta <- c(rep(FALSE, nNonRes), rep(TRUE, nRes))
    theta[(nNonRes + 1):(nNonRes + nRes)] <- 1
    trace_GLMM <- NULL
    for (i in 1:maxit) {
      deriv <- family$mu.eta(eta)
      z <- (eta - offset) + (y - mu)/deriv
      wGLM <- as.vector(deriv^2 / family$variance(mu))
      wGLM <- wGLM*w
      lRinv <- constructRinv(df = data, residual = residual, weights = wGLM)
      obj <- sparseMixedModels(y = z, X = Xs, Z = Z, lGinv = lGinv, lRinv = lRinv,
                               tolerance = tolerance, trace = trace, maxit = maxit,
                               theta = theta, fixedTheta = fixedTheta,
                               grpTheta = grpTheta)
      eta.old <- eta
      eta <- obj$yhat + offset
      mu <- family$linkinv(eta)
      theta <- obj$theta
      tol <- sum((eta - eta.old)^2) / sum(eta^2)

      aux.df <- data.frame(itOuter = rep(i, obj$nIter),
                           tol=c(rep(NA,obj$nIter-1), tol))
      trace_GLMM <- rbind(trace_GLMM, cbind(aux.df, obj$trace))
      if (trace) {
        cat("Generalized Linear Mixed Model iteration", i, ", tol=", tol, "\n")
      }
      if (tol < tolerance) {
        break;
      }
    }
    if (i == maxit) {
      warning("No convergence after ", maxit,
              " iterations of GLMM algorithm\n", call. = FALSE)
    }

    dev.residuals <- family$dev.resids(y, mu, w)
    deviance <- sum(dev.residuals)
  }
  ## Add names to ndx of coefficients.
  if (family$family == "multinomial") {
    ndxCf <- seq_along(1:n)
  } else {
    ndxCf <- seq_along(obj$a)
  }

  ## Fixed terms.
  ef <- cumsum(dim.f)
  sf <- ef - dim.f + 1
  ndxCoefF <- nameCoefs(coefs = ndxCf, desMat = X, termLabels = term.labels.f,
                        s = sf, e = ef, data = data, type = "fixed")
  ## Random terms.
  er <- sum(dim.f) + cumsum(dim.r)
  sr <- er - dim.r + 1
  ndxCoefR <- nameCoefs(coefs = ndxCf, termLabels = term.labels.r, s = sr,
                        e = er, data = data, group = group, type = "random")

  if (family$family == "multinomial") {
    dim.r <- sapply(dim.r, FUN= function(x) {x*nCat})
    dim.f <- sapply(dim.f, FUN= function(x) {x*nCat})

    #tmp1 <- rep(ndxCoefF, each=nCat)
    #tmp2 <- rep(respVar[-(nCat+1)], times=n)
    #ndxCoefF <- paste0(tmp1,"_",tmp2)

    #tmp1 <- rep(ndxCoefR, each=nCat)
    #tmp2 <- rep(respVar[-(nCat+1)], times=n)
    #ndxCoefR <- paste0(tmp1,"_",tmp2)
  }

  ## Combine result for fixed and random terms.
  ndxCoefTot <- c(ndxCoefF, ndxCoefR)
  names(ndxCoefTot) <- c(term.labels.f, term.labels.r)


  ## Nominal effective dimension for non-spline part
  if (nNonSplinesRandom > 0) {
    ## calculate NomEff dimension for non-spline part
    NomEffDimNonSplines <- calcNomEffDim(Xs, Z, dim.r[c(1:nNonSplinesRandom)], term.labels.r)
    ## combine with splines part
    NomEffDimRan <- c(NomEffDimNonSplines, NomEffDimRan)
  }
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
  if (family$family == "gaussian") {
    trace <- obj$trace
  } else {
    trace <- trace_GLMM
  }
  return(LMMsolveObject(logL = obj$logL,
                        sigma2e = obj$sigma2e,
                        tau2e = obj$tau2e,
                        EDdf = EDdf,
                        varPar = varPar,
                        VarDf = VarDf,
                        theta = obj$theta,
                        coefMME = obj$a,
                        ndxCoefficients = ndxCoefTot,
                        yhat = obj$yhat,
                        residuals = obj$residuals,
                        nIter = obj$nIter,
                        y = y,
                        X = Xs, # save spam format
                        Z = Z,
                        lGinv = lGinv,
                        lRinv = lRinv,
                        C = obj$C,
                        cholC = obj$cholC,
                        constantREML = constantREML,
                        dim = dim,
                        Nres = length(lRinv),
                        term.labels.f = term.labels.f,
                        term.labels.r = term.labels.r,
                        respVar = respVar,
                        splRes = splResList,
                        family = family,
                        deviance = deviance,
                        trace = trace))
}
