# Martin Boer, Biometris
#
rm(list=ls())
library(agridat)
library(asreml)
library(dplyr)
library(LMMsolver)
library(SpATS)
library(zoo)

# analysis of three models, within replicates:
# 1. Durban2003 data
# 2. Gilmour data, with row as fixed effect.
# 3. Gilmour data, with extra random term:

analysis = 3

if (analysis == 1) {
  data(durban.rowcol)
  dat <- durban.rowcol
  dat <- rename(dat, col=bed)
  dat <- mutate(dat, row = ((row-1) %% 8) + 1)
  fixed = as.formula(yield~rep+gen)
  extra_random_terms = FALSE
} else
{
  data(gilmour.serpentine)
  dat <- gilmour.serpentine
  dat <- mutate(dat, col = ((col-1) %% 5) + 1)

  if (analysis == 2) {
    # fixed effect for the model:
    fixed = as.formula(yield~rep+gen+R)
    extra_random_terms = FALSE
  } else {  # model 3
    fixed = as.formula(yield~rep+gen)
    extra_random_terms = TRUE
  }
}

dat <- transform(dat, R=factor(row), C=factor(col))
dat <- arrange(dat, rep, row, col)

head(dat)

devSAS <- function(obj.asr)
{
  Constant = log(2*pi)*summary(obj.asr)$nedf
  dev = -2*obj.asr$loglik + Constant
  dev
}

Nrep <- nlevels(dat$rep)
Nrow <- max(dat$row)
Ncol <- max(dat$col)

if (extra_random_terms == FALSE) {
  obj0.asr <- asreml(fixed, dat=dat)
} else {
  obj0.asr <- asreml(fixed, random=~rep:C+rep:R, dat=dat)
}
summary(obj0.asr)$varcomp
Constant = log(2*pi)*summary(obj0.asr)$nedf

devBaseline <- devSAS(obj0.asr)

if (extra_random_terms == FALSE)
{
  obj0.asr <- asreml(fixed, dat=dat, random=~units,
                   residual=~idv(rep):ar1(R):ar1(C),maxiter=100)
} else {
  obj0.asr <- asreml(fixed, dat=dat, random=~rep:C+rep:R+units,
                     residual=~idv(rep):ar1(R):ar1(C),maxiter=100)
}

summary(obj0.asr)$varcomp
devAR <- devSAS(obj0.asr)

# analysis LVxLV canonical form:

q1 <- Nrow
q2 <- Ncol

# first order, here we use scaling to stay consistent with results Emlyn:
Usc1 <- (1/sqrt(1/2))*calcUsc(q1, ord=1)
Usc2 <- (1/sqrt(1/2))*calcUsc(q2, ord=1)

lZ <- list()
lZ[[1]] = kronecker(diag(Nrep), kronecker(Usc1,rep(1,Ncol)))
lZ[[2]] = kronecker(diag(Nrep), kronecker(rep(1,Nrow),Usc2))
lZ[[3]] = kronecker(diag(Nrep), kronecker(Usc1, Usc2))

Z <- as.matrix(do.call("cbind", lZ))
colnames(Z) <- paste0("Z",1:ncol(Z))
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("row","col","rowcol"))

if (extra_random_terms == FALSE)
{
  obj1.LMM <- LMMsolve(fixed, randomMatrices=lM,data=dat_ext,eps=1.0e-8,
                     display=FALSE,monitor=FALSE)
} else {
  obj1.LMM <- LMMsolve(fixed, random=~rep:C+rep:R, randomMatrices=lM,data=dat_ext,eps=1.0e-8,
                       display=FALSE,monitor=FALSE)
}
obj1.LMM$ED
obj1.LMM$EDmax
obj1.LMM$sigma2e
obj1.LMM$tau2e

if (extra_random_terms == FALSE)  {
  obj1.asr <- asreml(fixed, random=~rep:grp(col)+rep:grp(row)+rep:grp(rowcol),
                   group=lM[], dat=dat_ext)
} else {
  obj1.asr <- asreml(fixed, random=~rep:C+rep:R+rep:grp(col)+rep:grp(row)+rep:grp(rowcol),
                     group=lM[], dat=dat_ext)
}
summary(obj1.asr)$varcomp
devSAS(obj1.asr)

obj1.LMM$logL - obj1.asr$loglik

invM_matrix <- function(n, type) {
  D = diff(diag(n), diff=1)
  Delta = 0.5*crossprod(D)
  if (type=='RWf') {
    alpha <- 1.0
    Con <- c(1,rep(0,n-1))
  } else if (type=='RWb') {
    alpha <- 1.0
    Con <- c(rep(0,n-1),1)
  } else if (type=="LV") {
    alpha <- 0.5
    Con <- c(sqrt(1/(n-1)),rep(0,n-2),sqrt(1/(n-1)))
  } else
  {
    stop("invM type not not defined")
  }
  M <- alpha*(2*Delta + tcrossprod(Con))
  rownames(M) <- paste0(1:n)
  colnames(M) <- colnames(1:n)
  attr(M, 'INVERSE') <- TRUE
  M
}

# define the random part of mixed model, and check the order
# for the interactions:
if (extra_random_terms == FALSE) {
  random = as.formula(~C:R:rep+C:rep+R:rep)
} else
{
  # to make distinction :
  dat$R2 <- dat$R
  dat$C2 <- dat$C
  random = as.formula(~C:R:rep + C:rep + R:rep + R2:rep + C2:rep)
}
mf <- model.frame(random, dat, drop.unused.levels = TRUE, na.action = NULL)
mt <- terms(mf)
attr(mt,"term.labels")

f.terms <- all.vars(mt)[attr(mt,"dataClasses") == "factor"]
Z1 <- model.matrix(mt, data = mf,
                   contrasts.arg = lapply(mf[,f.terms, drop = FALSE], contrasts, contrasts = FALSE))
head(colnames(Z1))

Mr_inv <- invM_matrix(Nrow,"LV")
Mc_inv <- invM_matrix(Ncol,"LV")
MrMc_inv = kronecker(Mr_inv, Mc_inv)

Mr_inv_ext <- kronecker.spam(diag.spam(Nrep), Mr_inv)
Mc_inv_ext <- kronecker.spam(diag.spam(Nrep), Mc_inv)
Mrc_inv_ext = kronecker.spam(diag.spam(Nrep), MrMc_inv)

# compare with LMMsolve using inverse matrix:
lGinv <- list()
lGinv[["C:rep"]] = as.spam(Mc_inv_ext)
lGinv[["R:rep"]] = as.spam(Mr_inv_ext)
lGinv[["C:R:rep"]] = as.spam(Mrc_inv_ext)
obj2.LMM <- LMMsolve(fixed, random=random, lGinverse=lGinv,data=dat, eps=1.0e-10,
                     display=TRUE)
obj2.LMM$logL
obj2.LMM$ED

# P-splines:
devLVxLVcano <- -2.0*obj1.LMM$logL + Constant
# LV x LV:
devLVxLV <- -2.0*obj2.LMM$logL + Constant

round(devBaseline,2)
round(devAR,2)
round(devLVxLV,2)
round(devLVxLVcano,2)



