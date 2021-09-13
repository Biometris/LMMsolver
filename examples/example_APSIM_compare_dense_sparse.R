#
# Martin Boer, Biometris, Wageningen, the Netherlands.
#
# example of analysis of APSIM data using splines.
#
rm(list = ls())
library(LMMsolver)
library(dplyr)
library(ggplot2)
library(zoo)
library(asreml)

# subset of APSIM simulation data set Daniela Bustos-Korts
dat_traj = read.csv('./examples/APSIM_Emerald.csv',header=TRUE, stringsAsFactors = FALSE)
sel_geno <- paste0("g",formatC(1:100, width=3,flag=0))
dat_traj = filter(dat_traj, year==2013, geno %in% sel_geno)

max_das <- max(dat_traj$das)

head(dat_traj)

# use ton's/ha
dat_traj$biomass <- dat_traj$biomass/1000

# data for the simulations:
xmin =  30
xmax = 130
step =   5
#step =   1
dat <- filter(dat_traj, das %in% seq(30, 130 , by=step))
dim(dat)

set.seed(1234)
dat$z <- dat$das
dat$env <- as.factor(dat$das)
N <- nrow(dat)
sigma2e = 0.1
#sigma2e = 0.00001
dat$ysim <- dat$biomass + rnorm(N,sd=sqrt(sigma2e))

labels <- unique(dat$geno)
dat$geno <- factor(dat$geno,levels=labels)
dat$env <- as.factor(dat$env)
dat$g_nr <- as.numeric(dat$geno)

Ngeno <- nlevels(dat$geno)
Nenv <- nlevels(dat$env)
Ngeno
Nenv

# here we define the splines
# 1: E cubical splines, second order differences
# 2: G first degree splines, ridge penalty

## 1) definitions for environmenal covariate:

degr1 = 3
pord1 = 2
xmin1 <- xmin
xmax1 <- xmax
nseg1 <- 15
#nseg1 <- 25

#dat$z <- as.numeric(scale(dat$z,scale=FALSE))

knots1 = PsplinesKnots(xmin1, xmax1, degr1, nseg1)
B1 <- Bsplines(knots1, dat$z)
q1 <- ncol(B1)

U1sc <- calcUsc(q1, pord1)

# definition of linear space, see
# ~gitprojects/MBnotes/sparse_mixed_model_splines.tex/pdf
nknots1 <- length(knots1)
tau <- rollmean(knots1[-c(1,nknots1)], k=degr1)
# B1(z) %*% tau = z
all.equal(as.vector(B1 %*% tau), dat$z)

## 1) definitions for genotype

# a very simple simple B-spline basis, helps that
# data can be in any order
knots2 = PsplinesKnots(1, Ngeno, 1, Ngeno-1)
B2 <- Bsplines(knots2, dat$g_nr)

# define matrix orthogonal to constant
J = matrix(data=1, nrow=Ngeno, ncol=Ngeno)
K = diag(Ngeno) - (1/Ngeno) * J
U2sc = eigen(K)$vectors[,-Ngeno, drop=FALSE]

# define the mixed model equations and solve:
lZ <- list()
lZ[[1]] = B1 %*% U1sc  #env
lZ[[2]] = B2 %*% U2sc  #geno
lZ[[3]] = RowKronecker(B1,B2) %*% kronecker(tau,  U2sc) # x.geno
lZ[[4]] = RowKronecker(B1,B2) %*% kronecker(U1sc, U2sc)
Z <- do.call("cbind", lZ)

dat_ext = cbind(dat, Z)
colnames(dat_ext) <- c(colnames(dat_ext)[1:ncol(dat)], paste0("col",1:ncol(Z)))

s <- proc.time()[3]
lM <- ndxMatrix(dat, lZ, c("f(z)","g","g.z","f_g(z)"))
obj.dense <- LMMsolve(ysim~z, group=lM, data=dat_ext,display=FALSE)
e <- proc.time()[3]
dense.time <- e-s
round(obj.dense$ED, 2)

lM <- ndxMatrix(dat, lZ, c("fz","g","gz","fgz"))
s <- proc.time()[3]
obj.asr <- asreml(ysim~z,random=~grp(fz)+grp(g)+grp(gz)+grp(fgz), group=lM[], data=dat_ext)
e <- proc.time()[3]
asr.time <- e-s
# sparse model:
D1 <- diff(diag(q1),    diff=2)
D2 <- diff(diag(Ngeno), diff=1)

lZ <- list()
lZ[[1]] = B1 %*% t(D1)  #env
lZ[[2]] = B2 %*% t(D2)  #geno
lZ[[3]] = RowKronecker(B1,B2) %*% kronecker(tau, t(D2))
lZ[[4]] = RowKronecker(B1,B2) %*% kronecker(t(D1), t(D2))

Z <- do.call("cbind", lZ)
dat_ext = cbind(dat, Z)

lM <- ndxMatrix(dat, lZ, c("f(z)","g","g.z","f_g(z)"))
I_g <- diag(1,Ngeno)
DtD1 <- crossprod(D1)
lGinv <- list()
lGinv[['f(z)']] <- as.spam(D1 %*% DtD1 %*% t(D1))
lGinv[['g']]   <- as.spam(D2 %*% I_g %*% t(D2))
lGinv[['g.z']] <- as.spam(D2 %*% I_g %*% t(D2))
lGinv[['f_g(z)']] <- as.spam(kronecker(D1 %*% DtD1 %*% t(D1), D2 %*% I_g %*% t(D2)))
names(lGinv)
s <- proc.time()[3]
obj.sparse <- LMMsolve(ysim~z, group=lM,lGinverse=lGinv, data=dat_ext,
                         display=TRUE,trace=FALSE)
e <- proc.time()[3]
sparse.time <- e-s
round(obj.sparse$ED, 2)

# compare logL:
obj.asr$loglik
obj.dense$logL
obj.sparse$logL

# compare ED's
obj.dense$ED
obj.sparse$ED

# compare some coefficients:
fz1 <- U1sc %*% coef(obj.dense)$'f(z)'
fz2 <- t(D1) %*% coef(obj.sparse)$'f(z)'
all.equal(fz1, fz2)

fgz1 <- kronecker(U1sc, U2sc) %*% coef(obj.dense)$'f_g(z)'
fgz2 <- kronecker(t(D1),t(D2)) %*% coef(obj.sparse)$'f_g(z)'
all.equal(fgz1, fgz2)

dense.time
sparse.time
asr.time
