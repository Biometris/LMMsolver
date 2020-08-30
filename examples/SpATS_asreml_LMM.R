library(agridat)
library(asreml)
library(dplyr)
library(LMMsolver)
library(SpATS)
library(zoo)

devSAS <- function(obj.asr)
{
  Constant = log(2*pi)*summary(obj.asr)$nedf
  dev = -2*obj.asr$loglik + Constant
  dev
}

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

obj0.SpATS <- SpATS(response="yield",spatial=~PSANOVA(col,row, nseg=c(Ncols-1,Nrows-1)),
                   genotype="gen",fixed=~rep, data=dat, control=list(monitor=0,tolerance=1.0e-6))
summary(obj0.SpATS)
plot(obj0.SpATS)

# model with row, col corrections:

obj1.SpATS <- SpATS(response="yield",spatial=~PSANOVA(col, row, nseg=c(Ncols-1,Nrows-1)),
                   genotype="gen",fixed=~rep, random=~R+C, data=dat,
                   control=list(monitor=0))
summary(obj1.SpATS)


# Ok, same result as Hans-Peter:
obj0.asr <- asreml(yield~rep+gen, data=dat,trace=FALSE)
devSAS(obj0.asr)

# Ok, same result as Hans-Peter:
obj1.asr <- asreml(yield~rep+gen, data=dat,random=~units, resid = ~id(rep):ar1(Rw):ar1(C),trace=FALSE)
summary(obj1.asr)$varcomp
devSAS(obj1.asr)

# Ok, same result as Hans-Peter:
obj2.asr <- asreml(yield~rep+gen, data=dat,random=~units, resid = ~ar1(R):ar1(C),trace=FALSE)
summary(obj2.asr)$varcomp
devSAS(obj2.asr)

# definition of splines...

x1_sc <- scale(dat$row)
x2_sc <- scale(dat$col)

knots1 <- PsplinesKnots(1,Nrows, degree=3, nseg=Nrows-1)
knots2 <- PsplinesKnots(1,Ncols, degree=3, nseg=Ncols-1)

tau1 <- rollmean(knots1[-c(1,length(knots1))], k=3)
tau2 <- rollmean(knots2[-c(1,length(knots2))], k=3)

B1 <- Bsplines(knots1, dat$row)
B2 <- Bsplines(knots2, dat$col)

# B1(x) tau1 = x, B2(x) tau2 = x:
range(B1 %*% tau1 - dat$row)
range(B2 %*% tau2 - dat$col)

q1 <- ncol(B1)
q2 <- ncol(B2)

# first order:
Usc1 <- calcUsc(q1, ord=1)
Usc2 <- calcUsc(q2, ord=1)

lZ <- list()
lZ[[1]] = B1 %*% Usc1  
lZ[[2]] = B2 %*% Usc2 
lZ[[3]] = RowKronecker(B1, B2) %*% kronecker(Usc1, Usc2)

Z <- as.matrix(do.call("cbind", lZ))
colnames(Z) <- paste0("Z",1:ncol(Z))
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("row","col","rowcol"))

obj1.LMM <- LMMsolve(yield~rep+gen, randomMatrices=lM,data=dat_ext,eps=1.0e-4,
                display=TRUE,monitor=TRUE)
obj1.LMM$ED
obj1.LMM$EDmax

obj3.asr <- asreml(yield~rep+gen, random=~grp(row)+grp(col)+grp(rowcol),
                   group=lM[], dat=dat_ext)
summary(obj3.asr)$varcomp

devSAS(obj3.asr)

# second order:
Usc1 <- calcUsc(q1, ord=2)
Usc2 <- calcUsc(q2, ord=2)

lZ <- list()
lZ[[1]] = B1 %*% Usc1  
lZ[[2]] = B2 %*% Usc2 
lZ[[3]] = RowKronecker(B1, B2) %*% kronecker(Usc1, tau2)
lZ[[4]] = RowKronecker(B1, B2) %*% kronecker(tau1, Usc2)
lZ[[5]] = RowKronecker(B1, B2) %*% kronecker(Usc1, Usc2)

Z <- as.matrix(do.call("cbind", lZ))
colnames(Z) <- paste0("Z",1:ncol(Z))
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("fr","fc","fr.c","r.fc","fr.fc"))

obj2.LMM <- LMMsolve(yield~rep+gen+row+col+row:col, randomMatrices=lM,data=dat_ext,eps=1.0e-6,
                     display=TRUE,monitor=TRUE)
obj2.LMM$ED
obj2.LMM$EDmax

obj4.asr <- asreml(yield~rep+gen+row+col+row:col, 
                   random=~grp(fr)+grp(fc)+grp(fr.c)+grp(r.fc)+grp(fr.fc),
                   group=lM[], dat=dat_ext)
obj4.asr <- update(obj4.asr)
summary(obj4.asr)$varcomp

Constant_Psplines2order = log(2*pi)*summary(obj4.asr)$nedf

devSAS(obj4.asr)

obj2.LMM$logL - obj4.asr$loglik

# second order, other model:
#Usc1 <- calcUsc(q1, ord=2)
#Usc2 <- calcUsc(q2, ord=2)

ord=2
D1 <- diff(diag(q1), diff=2)
DtD1 <- crossprod(D1)
U1 <- eigen(DtD1)$vectors[,1:(q1-ord)]
d1 <- eigen(DtD1)$values[1:(q1-ord)]

D2 <- diff(diag(q2), diff=2)
DtD2 <- crossprod(D2)
U2 <- eigen(DtD2)$vectors[,1:(q2-ord)]
d2 <- eigen(DtD2)$values[1:(q2-ord)]

X1 <- obj0.SpATS$terms$spatial$MM$MM1$X
X2 <- obj0.SpATS$terms$spatial$MM$MM2$X
X12 <- RowKronecker(X1,X2)
dat$x1 <- X12[,2]
dat$x2 <- X12[,3]
dat$x12 <- X12[,4]

lZ <- list()
lZ[[1]] = B1 %*% Usc1  
lZ[[2]] = B2 %*% Usc2 
#lZ[[3]] = RowKronecker(B1, B2) %*% kronecker(Usc1, tau2)
#lZ[[4]] = RowKronecker(B1, B2) %*% kronecker(tau1, Usc2)
lZ[[3]] = (B1 %*% Usc1) * dat$x2
lZ[[4]] = dat$x1 * (B2 %*% Usc2)

pord = 2
d3 <- c(d1 %x% rep(1, q2 - pord) + rep(1, q1 - pord) %x% d2)
#d3 <- c(rep(1, c2 - pord) %x% d1 + d2 %x% rep(1, c1 - pord))

# here we use U instead of Usc:
lZ[[5]] = RowKronecker(B1, B2) %*% kronecker(U1,U2) %*% diag(1/sqrt(d3))


Z <- as.matrix(do.call("cbind", lZ))
colnames(Z) <- paste0("Z",1:ncol(Z))
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("fr","fc","fr.c","r.fc","fr.fc"))

lGinv <- list()
#lGinv[['fr.fc']] <- kronecker.spam(diag.spam(d1),diag.spam(d2))

# see SpATS paper: 
#lGinv[['fr.fc']] <- kronecker.spam(diag.spam(d1),diag.spam(q2-2)) +
#                    kronecker.spam(diag.spam(q1-2),diag.spam(d2))
    
obj3.LMM <- LMMsolve(yield~rep+gen+x1+x2+x12, randomMatrices=lM,
                     data=dat_ext,eps=1.0e-6,
                     display=TRUE,monitor=TRUE)

obj2.LMM$logL
obj3.LMM$logL

obj5.asr <- asreml(yield~rep+gen+x1+x2+x12, 
                   random=~grp(fr)+grp(fc)+grp(fr.c)+grp(r.fc)+grp(fr.fc),
                   group=lM[], dat=dat_ext)

devSAS(obj5.asr)
dev = -2*obj3.LMM$logL + Constant_Psplines2order
dev

logL.SpATS <- -0.5*obj0.SpATS$deviance
obj3.LMM$logL
logL.SpATS
obj5.asr$loglik

logL.SpATS-obj3.LMM$logL

