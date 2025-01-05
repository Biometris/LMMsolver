fitLMM <- function(y, X, Z, w, lGinv, tolerance, trace, maxit,
              theta, grpTheta, family, offset, dim.f, dim.r,
              term.labels.f, term.labels.r, respVar, NomEffDimRan,
              varPar, splResList, residual, group,
              nNonSplinesRandom, scFactor, data) {
  ## Convert to spam matrix and cleanup
  Xs <- spam::as.spam.dgCMatrix(X)
  Xs <- spam::cleanup(Xs)
  npar <- 0
  if (!is.null(dim.f)) {
    npar <- npar + sum(dim.f)
  }
  if (!is.null(dim.r)) {
    npar <- npar + sum(dim.r)
  }


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
  } else if (family$family != "multinomial") {
    ## MB, 23 jan 202
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
  } else {
    ## multinomial family
    YY <- y
    n <- nrow(YY)
    nCat <- length(respVar) - 1
    sY <- rowSums(YY)
    YY <- YY/sY

    ## allow for incorrect scores, to make the algorithm more stable
    YY <- t(apply(YY, MARGIN=1, FUN=function(x) {
      eps <- 1.0e-6
      ndx_zero <- which(x==0)
      ndx_pos  <- which(x>0)
      x[ndx_zero] <- x[ndx_zero] + eps
      nZeros <- length(ndx_zero)
      nPos <- length(ndx_pos)
      x[ndx_pos] <- x[ndx_pos] - (nZeros/nPos)*eps
      return(x)
    }))

    YY <- YY[,-(nCat+1)] # remove last column
    y <- as.vector(t(YY))
    Xs <- Xs %x% spam::diag.spam(nCat)
    Z <- Z %x% spam::diag.spam(nCat)
    w <- rep(w, nCat)
    lGinv <- lapply(lGinv, FUN = function(x) {x %x% spam::diag.spam(nCat)})

    # remark: We can add some prior information to initial
    # values based on the scores:
    eta <- rep(0, n*nCat)

    # Initialize the sparse block structure for W and Dinv
    W <- spam::diag.spam(1, n) %x% spam::spam(x=1, nrow=nCat, ncol=nCat)
    Dinv <- spam::diag.spam(1, n) %x% spam::spam(x=1, nrow=nCat, ncol=nCat)

    #theta <- c(1,1)                # check this!
    #fixedTheta <- c(FALSE, TRUE)   # check this!
    nRes <- length(lRinv)
    nNonRes <- length(theta) - nRes
    fixedTheta <- c(rep(FALSE, nNonRes), rep(TRUE, nRes))
    theta[(nNonRes + 1):(nNonRes + nRes)] <- 1

    trace_GLMM <- NULL
    for (i in 1:maxit) {
      Eta <- matrix(data=eta, ncol=nCat, nrow=n, byrow=TRUE)
      Pi <- t(apply(X=Eta, MARGIN=1, FUN = family$linkinv))
      pi_vec <- as.vector(t(Pi))

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

      z <- eta + Dinv %*% (y - pi_vec)
      lRinv <- list(residual = W)
      attr(lRinv, "cnt") <- n*nCat # correct?
      obj <- sparseMixedModels(z, X = Xs, Z = Z,
                               lGinv = lGinv, lRinv = lRinv, trace=trace,
                               fixedTheta = fixedTheta, theta = theta,
                               tolerance = tolerance)
      eta_old  <- eta
      eta <- obj$yhat
      theta <- obj$theta
      tol <- sum((eta - eta_old)^2)/sum(eta^2)

      aux.df <- data.frame(itOuter = rep(i, obj$nIter),
                           tol=c(rep(NA,obj$nIter-1), tol))
      trace_GLMM <- rbind(trace_GLMM, cbind(aux.df, obj$trace))
      if (trace) {
        cat("Generalized Linear Mixed Model iteration", i, ", tol=", tol, "\n")
      }
      #cat("iter ", i, "     ", tol, "    ", obj$logL, "   ", obj$ED, "\n")
      if (tol < 1.0e-10) {
        if (trace) {
          cat("convergence after", i, "iterations\n")
        }
        break;
      }
    }
    mu <- pi_vec
    dev.residuals <- family$dev.resids(y, mu, w)
    deviance <- sum(dev.residuals)
  }
  ndxCf <- seq_len(npar)

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
  }

  ## Combine result for fixed and random terms.
  ndxCoefTot <- c(ndxCoefF, ndxCoefR)
  names(ndxCoefTot) <- c(term.labels.f, term.labels.r)

  ## for multinomial extra coefficients and names
  if (family$family== "multinomial") {
    ndxCoefTot <- extend_coef(ndxCoefTot, respVar)
  }

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
  N <- nrow(Xs)
  p <- ncol(Xs)
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
