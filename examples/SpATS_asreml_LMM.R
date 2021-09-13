# Martin Boer, Biometris
rm(list=ls())
library(agridat)
library(asreml)
library(dplyr)
library(LMMsolver)
library(SpATS)
library(zoo)

scale_fixed_effect = TRUE

# function to calculate deviance:
devSAS <- function(obj.asr)
{
  Constant = log(2*pi)*summary(obj.asr)$nedf
  dev = -2*obj.asr$loglik + Constant
  dev
}

data(durban.rowcol)

dat <- durban.rowcol
dat <- rename(dat, col=bed)

dat <- dat[order(dat$row, dat$col),]
head(dat)

Nrows <- max(dat$row)
Ncols <- max(dat$col)

nsegR <- Nrows/2
nsegC <- Ncols/2
#nsegR <- Nrows-1
#nsegC <- Ncols-1
obj0.SpATS <- SpATS(response="yield",spatial=~PSANOVA(col,row, nseg=c(nsegC,nsegR)),
                   genotype="gen",fixed=~rep, data=dat, control=list(trace=0,tolerance=1.0e-10))
summary(obj0.SpATS)
plot(obj0.SpATS)

# obtain fixed part from SpATS:
if (scale_fixed_effect) {
  X1 <- obj0.SpATS$terms$spatial$MM$MM1$X
  X2 <- obj0.SpATS$terms$spatial$MM$MM2$X
  X12 <- RowKronecker(X1,X2)
  dat$x1 <- X12[,2]
  dat$x2 <- X12[,3]
  dat$x12 <- X12[,4]
  #dat$x1 <- scale(dat$row)
  #dat$x2 <- scale(dat$col)
  #dat$x12 <- dat$x1 * dat$x2
} else {
dat$x1 <- dat$row
dat$x2 <- dat$col
dat$x12 <- dat$x1 * dat$x2
}

tolerance = 1.0e-14
xmin1 = min(dat$x1) - tolerance
xmax1 = max(dat$x1) + tolerance
xmin2 = min(dat$x2) - tolerance
xmax2 = max(dat$x2) + tolerance

# definition of splines...
knots1 <- PsplinesKnots(xmin1,xmax1, degree=3, nseg=nsegR)
knots2 <- PsplinesKnots(xmin2,xmax2, degree=3, nseg=nsegC)

B1 <- Bsplines(knots1, dat$x1)
B2 <- Bsplines(knots2, dat$x2)

q1 <- ncol(B1)
q2 <- ncol(B2)

one1 <- matrix(data=1, ncol=1, nrow=q1)
one2 <- matrix(data=1, ncol=1, nrow=q2)

range(B1 %*% one1)
range(B2 %*% one2)

tau1 <- rollmean(knots1[-c(1,length(knots1))], k=3)
tau2 <- rollmean(knots2[-c(1,length(knots2))], k=3)

# B1(x1) tau1 = x1, B2(x2) tau2 = x2:
range(B1 %*% tau1 - dat$x1)
range(B2 %*% tau2 - dat$x2)

# second order:
pord = 2
D1 <- diff(diag(q1), diff=pord)
DtD1 <- crossprod(D1)
U1 <- eigen(DtD1)$vectors[,1:(q1-pord)]
d1 <- eigen(DtD1)$values[1:(q1-pord)]

D2 <- diff(diag(q2), diff=pord)
DtD2 <- crossprod(D2)
U2 <- eigen(DtD2)$vectors[,1:(q2-pord)]
d2 <- eigen(DtD2)$values[1:(q2-pord)]

# see SpATS paper, used for interaction term:
d3 <- c(d1 %x% rep(1, q2 - pord) + rep(1, q1 - pord) %x% d2)
objb.LMM <- LMMsolve(yield~rep+gen+x1+x2+x12,
                     data=dat,tolerance=1.0e-10,
                     display=FALSE,trace=FALSE)

objb.asr <- asreml(yield~rep+gen+x1+x2+x12, dat=dat)

logL.SpATS <- -0.5*obj0.SpATS$deviance

devSAS(objb.asr)


lZ <- list()
# we have to calculate RowKronecker product only once:
B12 <- RowKronecker(B1, B2)
lZ[[1]] = B12 %*% kronecker(U1, one2) %*% diag(1/sqrt(d1))
lZ[[2]] = B12 %*% kronecker(one1, U2) %*% diag(1/sqrt(d2))
lZ[[3]] = B12 %*% kronecker(U1, tau2) %*% diag(1/sqrt(d1))
lZ[[4]] = B12 %*% kronecker(tau1, U2) %*% diag(1/sqrt(d2))
lZ[[5]] = B12 %*% kronecker(U1,   U2) %*% diag(1/sqrt(d3))

Z <- as.matrix(do.call("cbind", lZ))
colnames(Z) <- paste0("Z",1:ncol(Z))
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("fr","fc","fr.c","r.fc","fr.fc"))

obj0.LMM <- LMMsolve(yield~rep+gen+x1+x2+x12, group=lM,
                     data=dat_ext,tolerance=1.0e-10,
                     display=FALSE,trace=FALSE)

obj0.asr <- asreml(yield~rep+gen+x1+x2+x12,
                   random=~grp(fr)+grp(fc)+grp(fr.c)+grp(r.fc)+grp(fr.fc),
                   group=lM[], dat=dat_ext)

logL.SpATS <- -0.5*obj0.SpATS$deviance

logL.SpATS
obj0.LMM$logL
obj0.asr$loglik

obj0.LMM$logL - logL.SpATS
obj0.asr$loglik - logL.SpATS

devSAS(obj0.asr)

# model G6 (Wood formulation):
# first order:
Usc1 <- calcUsc(q1, ord=2)
Usc2 <- calcUsc(q2, ord=2)

lZ <- list()
# we have to calculate RowKronecker product only once:
B12 <- RowKronecker(B1, B2)
lZ[[1]] = B12 %*% kronecker(Usc1, one2)
lZ[[2]] = B12 %*% kronecker(one1, Usc2)
lZ[[3]] = B12 %*% kronecker(Usc1, tau2)
lZ[[4]] = B12 %*% kronecker(tau1, Usc2)
lZ[[5]] = B12 %*% kronecker(Usc1, Usc2)

range(t(kronecker(Usc1, one2)) %*% kronecker(Usc1,tau2))

Z <- as.matrix(do.call("cbind", lZ))
colnames(Z) <- paste0("Z",1:ncol(Z))
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("fr","fc","fr.c","r.fc","fr.fc"))

obj1.LMM <- LMMsolve(yield~rep+gen+x1+x2+x12, group=lM,
                     data=dat_ext,tolerance=1.0e-10,
                     display=FALSE,trace=FALSE)

obj1.asr <- asreml(yield~rep+gen+x1+x2+x12,
                   random=~grp(fr)+grp(fc)+grp(fr.c)+grp(r.fc)+grp(fr.fc),
                   group=lM[], dat=dat_ext)

obj1.LMM$logL
obj1.asr$loglik

devSAS(obj1.asr)


# Model without interaction

lZ <- list()
# we have to calculate RowKronecker product only once:
B12 <- RowKronecker(B1, B2)
lZ[[1]] = B12 %*% kronecker(Usc1, one2)
lZ[[2]] = B12 %*% kronecker(one1, Usc2)
lZ[[3]] = B12 %*% kronecker(Usc1, tau2)
lZ[[4]] = B12 %*% kronecker(tau1, Usc2)

Z <- as.matrix(do.call("cbind", lZ))
colnames(Z) <- paste0("Z",1:ncol(Z))
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("fr","fc","fr.c","r.fc"))

obj2.LMM <- LMMsolve(yield~rep+gen+x1+x2+x12, group=lM,
                     data=dat_ext,tolerance=1.0e-10,
                     display=FALSE,trace=FALSE)

obj2.asr <- asreml(yield~rep+gen+x1+x2+x12,
                   random=~grp(fr)+grp(fc)+grp(fr.c)+grp(r.fc),
                   group=lM[], dat=dat_ext)

obj2.LMM$logL
obj2.asr$loglik


devSAS(objb.asr) # base model
devSAS(obj2.asr) # without G5 or G6 interaction
devSAS(obj0.asr)  # SpATS, G5 model
devSAS(obj1.asr)  # G6 interaction

Constant = log(2*pi)*summary(obj1.asr)$nedf
dev = -2*obj1.LMM$logL + Constant
dev

devSpATS = obj0.SpATS$deviance + Constant
devSpATS
