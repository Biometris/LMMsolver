library(asreml)
library(dplyr)
library(LMMsolver)
data(john.alpha, package = "agridat")
head(john.alpha)
dat <- john.alpha
Nblk <- nlevels(dat$block)*nlevels(dat$rep)
dat$block2 <- as.factor(rep(paste0("BLK",  formatC(1:Nblk,width=2,flag="0")),each=4))

obj0 <- asreml(fixed = yield ~ gen+rep,
               random = ~block:at(rep),
               data = dat, tolerance=1.0e-10)
summary(obj0)

obj1 <- asreml(fixed = yield ~ gen+rep,
               random = ~at(rep,"R1"):block + at(rep,"R2"):block,
               data = dat, tolerance=1.0e-10)
summary(obj1)

# by foot in LMMsolver
dat_ext <- dat
dat_ext$R1_cond <- 1*(dat$rep=="R1")
dat_ext$R2_cond <- 1*(dat$rep=="R2")

obj2 <- LMMsolve(fixed=yield~gen+rep,
                random=~block:R1_cond + block:R2_cond,data=dat_ext)
summary(obj2)
obj1$loglik
obj2$logL

# define cf() function, to allow for variances conditional on factor:
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
  ndx <- which(colSums(Z)!=0)
  Z <- Z[, ndx]
  Z <- spam::as.spam.dgCMatrix(Z)
  termlabel <- paste0("cf(",cName,", ",level,"):",vName)
  cflabel <- paste(termlabel, levels(var)[ndx],sep="_")
  L <- list(Z=Z, termlabel = termlabel, cflabel=cflabel, dim.r=ncol(Z))
  return(L)
}

CondFactor <- function(random, data) {
  tf <- terms.formula(random, specials = c("cf"))
  f <- attr(tf, "term.labels")
  ndxAt <- attr(tf,"specials")$cf
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
  return(list(Z=Z, termlabel=termlabel, cflabel=cflabel, dim.r=dim.r, random=random))
}

random = ~block + cf(block2,rep,"R1") + cf(block2,rep,"R2")

obj <- CondFactor(random, dat)
dim(obj$Z)
obj$dim.r

