# Martin Boer, Biometris
rm(list=ls())
library(agridat)
library(asreml)
library(dplyr)
library(LMMsolver)
library(SpATS)
library(zoo)

data(durban.rowcol)

dat <- durban.rowcol
dat <- rename(dat, col=bed)
dat <- transform(dat, rowW = ((row-1) %% 8) + 1)
head(dat)

dat <- transform(dat, R=factor(row), C=factor(col), Rw=factor(rowW))
dat <- dat[order(dat$R, dat$C),]
head(dat)

Nrows <- max(dat$row)
Ncols <- max(dat$col)
# model without row/col corrections:

nsegR <- Nrows/2
nsegC <- Ncols/2
obj0.SpATS <- SpATS(response="yield",spatial=~PSANOVA(col,row, nseg=c(nsegC,nsegR)),
                   genotype="gen",fixed=~rep, data=dat, control=list(monitor=0,tolerance=1.0e-6))
summary(obj0.SpATS)
plot(obj0.SpATS)

# obtain fixed part from SpATS:
X1 <- obj0.SpATS$terms$spatial$MM$MM1$X
X2 <- obj0.SpATS$terms$spatial$MM$MM2$X
X12 <- RowKronecker(X1,X2)
dat$x1 <- X12[,2]
dat$x2 <- X12[,3]
dat$x12 <- X12[,4]

# definition of splines...
knots1 <- PsplinesKnots(1,Nrows, degree=3, nseg=nsegR)
knots2 <- PsplinesKnots(1,Ncols, degree=3, nseg=nsegC)

tau1 <- rollmean(knots1[-c(1,length(knots1))], k=3)
tau2 <- rollmean(knots2[-c(1,length(knots2))], k=3)

B1 <- Bsplines(knots1, dat$row)
B2 <- Bsplines(knots2, dat$col)

# B1(x) tau1 = x, B2(x) tau2 = x:
range(B1 %*% tau1 - dat$row)
range(B2 %*% tau2 - dat$col)

q1 <- ncol(B1)
q2 <- ncol(B2)

# second order:
pord = 2
Usc1 <- calcUsc(q1, ord=pord)
Usc2 <- calcUsc(q2, ord=pord)

D1 <- diff(diag(q1), diff=2)
DtD1 <- crossprod(D1)
U1 <- eigen(DtD1)$vectors[,1:(q1-pord)]
d1 <- eigen(DtD1)$values[1:(q1-pord)]

D2 <- diff(diag(q2), diff=2)
DtD2 <- crossprod(D2)
U2 <- eigen(DtD2)$vectors[,1:(q2-pord)]
d2 <- eigen(DtD2)$values[1:(q2-pord)]

# see SpATS paper, used for interaction term:
d3 <- c(d1 %x% rep(1, q2 - pord) + rep(1, q1 - pord) %x% d2)

lZ <- list()
lZ[[1]] = B1 %*% Usc1
lZ[[2]] = B2 %*% Usc2
#lZ[[3]] = RowKronecker(B1, B2) %*% kronecker(Usc1, tau2)
#lZ[[4]] = RowKronecker(B1, B2) %*% kronecker(tau1, Usc2)
lZ[[3]] = (B1 %*% Usc1) * dat$x2
lZ[[4]] = dat$x1 * (B2 %*% Usc2)
lZ[[5]] = RowKronecker(B1, B2) %*% kronecker(U1,U2) %*% diag(1/sqrt(d3))

Z <- as.matrix(do.call("cbind", lZ))
colnames(Z) <- paste0("Z",1:ncol(Z))
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("fr","fc","fr.c","r.fc","fr.fc"))

obj0.LMM <- LMMsolve(yield~rep+gen+x1+x2+x12, randomMatrices=lM,
                     data=dat_ext,eps=1.0e-6,
                     display=TRUE,monitor=TRUE)

obj0.asr <- asreml(yield~rep+gen+x1+x2+x12,
                   random=~grp(fr)+grp(fc)+grp(fr.c)+grp(r.fc)+grp(fr.fc),
                   group=lM[], dat=dat_ext)

logL.SpATS <- -0.5*obj0.SpATS$deviance

logL.SpATS
obj0.LMM$logL
obj0.asr$loglik

obj0.LMM$logL - logL.SpATS
obj0.asr$loglik - logL.SpATS
