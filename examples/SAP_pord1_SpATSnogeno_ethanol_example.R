#' ---
#' title: SAP example sparse
#' author: Martin Boer, Biometris, WUR, Wageningen, the Netherlands
#' ---


rm(list=ls())

library(spam)
library(dplyr)
library(SpATS)
library(LMMsolver)
library(ggplot2)
library(JOPS)

# Get the data
library(SemiPar)
data(ethanol)

#' SAP model using SpATS.nogeno
#' ==================

nseg <- c(25, 25)
degr <- 3
pord <- 1

# Fit the PS-ANOVA model
obj1 <- SpATS.nogeno(response = "NOx",
                     spatial = ~SAP(E, C,  nseg = nseg, pord=pord),
                     data = ethanol,
                     control = list(maxit = 100, tolerance = 1.0e-10,
                                    monitoring = 2, update.psi = FALSE))

summary(obj1)

# Compute component surface and their sum on a fine grid
Tr = obtain.spatialtrend(obj1, grid = c(100, 100))

# Plot surface and contours
image(Tr$row.p, Tr$col.p, Tr$fit, col = terrain.colors(100), xlab = 'C', ylab = 'E')
contour(Tr$row.p, Tr$col.p, Tr$fit, add = TRUE, col = 'blue')
points(ethanol$C, ethanol$E, pch = '+')

#' solve using sparse formulation SAP
#' ==================

dat <- ethanol
x1 <- dat$E
x2 <- dat$C

knots1 <- PsplinesKnots(min(x1), max(x1), degree=degr, nseg=nseg[1])
knots2 <- PsplinesKnots(min(x2), max(x2), degree=degr, nseg=nseg[2])

B1 <- Bsplines(knots1, x1)
B2 <- Bsplines(knots2, x2)
q1 <- ncol(B1)
q2 <- ncol(B2)

D1 <- diff(diag(q1), diff=pord)
DtD1 <- crossprod(D1)

D2 <- diff(diag(q2), diff=pord)
DtD2 <- crossprod(D2)

# we have to calculate RowKronecker product only once:
B21 <- RowKronecker(B2, B1)

Nind <- nrow(dat)
Nblocks <- nlevels(dat$plotId)
X<-matrix(data=1,ncol=1,nrow=Nind)

# make Rinv
#Nind <- nrow(X)
lRinv = list()
lRinv[[1]] = diag.spam(1,Nind)
names(lRinv) = "residual"

Z <- B21

dim.r = c(ncol(B21), Nblocks)
e <- cumsum(dim.r)
s <- e - dim.r + 1

# change the constraints, rest same as before..
C2 <- c(1, rep(0,q1*q2-2), 1)
CCt2 <- tcrossprod(C2)
lGinv <- list()
lGinv[[1]] <- bdiag.spam(as.spam(kronecker.spam(diag(q2), DtD1) + CCt2))
lGinv[[2]] <- bdiag.spam(as.spam(kronecker.spam(DtD2, diag(q1)) + CCt2))

y <- dat$NOx
s <- proc.time()[3]
obj2 = sparseMixedModels(y,X,Z,lGinv,lRinv,maxiter=300,
                          eps=1.0e-10,display=TRUE,monitor=TRUE)
e <- proc.time()[3]
cat("Computation time:", e-s, "seconds")
obj2$ED

p = ncol(X)
# sum of first and last element equal to zero, as for 1D LV model:
obj2$a[p+1]
obj2$a[p+ncol(B21)]

#' compare results
#' ==================

# compare SAP SpATS with LMMsolver..
dev1 <- as.numeric(obj1$deviance)
dev2 <- -2*obj2$logL
dev1
dev2

dev1-dev2

