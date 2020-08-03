#' ---
#' title: Spatial analysis using first order P-splines (Whittaker)
#' author: Martin Boer, Biometris, WUR
#' ---

#' Load Gilmour 1997 data
#' ==================

rm(list=ls())

library(asreml)
library(dplyr)
library(fields)
library(LMMsolver)

df = read.csv("Piepho_Williams2010_example2.csv",na.string='.')

# 5 columns within replicate
nrows = 22; ncols=5;
nrep = 3
ngeno = 107

# 3 genotypes twice in each replicate, see Gilmour 1997
nrep*(ngeno-3) + 2*nrep*3

nobs <- nrow(df)

df$Row = as.factor(df$row)
df$Col = as.factor(df$col)
df$gen = as.factor(df$gen)
df$rep = as.factor(df$rep)

## baseline model using asreml, row and column defined within replicate:
obj1 = asreml(fixed = yield ~ gen + rep, random = ~rep:Row+rep:Col, data = df, trace = FALSE)
Constant = log(2*pi)*(obj1$nedf)
obj1$loglik

## same model, using LMMsolve:
obj2 = LMMsolve(fixed=yield~-1+gen+rep, random=~rep:Row+rep:Col, data=df)
obj2$logL

#' P-splines model, baseline + P(row)
#' ==================

U1sc <- calcUsc(nrows, ord=1)
tmp1 <- kronecker(U1sc, rep(1,ncols))
lZ <- list()
lZ[[1]] <- kronecker(diag(nrep),tmp1)
Z <- do.call("cbind",lZ)
df_ext = cbind(df, lZ)
lM <- ndxMatrix(df, lZ, "Prow")
obj3 = LMMsolve(fixed=yield~-1+gen+rep, random=~rep:Row+rep:Col,randomMatrices=lM, data=df_ext)
obj3$logL
obj3$ED

# sum equal to zero, follows by definition from U1sc:
A1 <- kronecker(diag(nrep), rep(1,nrows))
B1 <- kronecker(diag(nrep), U1sc)
C1 <- t(A1) %*% B1
range(C1)

C1 %*% coef(obj3)$Prow

#' P-splines model, baseline + P(row) + P(col)
#' ==================

U2sc <- calcUsc(ncols, ord=1)
tmp2 <- kronecker(rep(1,nrows), U2sc)
lZ[[2]] = kronecker(diag(nrep),tmp2)

Z <- do.call("cbind",lZ)
df_ext = cbind(df, lZ)
lM <- ndxMatrix(df, lZ, c("Prow","Pcol"))
obj4 = LMMsolve(fixed=yield~-1+gen+rep, random=~rep:Row+rep:Col,randomMatrices=lM, data=df_ext)
obj4$logL
obj4$ED

dev4 = -2*obj4$logL + Constant
dev4

# sum equal to zero, follows by definition from U1sc
C1 %*% coef(obj4)$Prow

# sum equal to zero, follows by definition from U2sc:
A2 <- kronecker(diag(nrep), rep(1, ncols))
B2 <- kronecker(diag(nrep), U2sc)
C2 <- t(A2) %*% B2
range(C2)
C2 %*% coef(obj4)$Pcol

#' P-splines model, baseline + P(row) + P(col) + P(row):P(col)
#' ==================

U12sc = kronecker(U1sc, U2sc)
lZ[[3]] = kronecker(diag(nrep), U12sc)

Z <- do.call("cbind",lZ)
df_ext = cbind(df, lZ)
lM <- ndxMatrix(df, lZ, c("Prow","Pcol","ProwxPcol"))
obj5 = LMMsolve(fixed=yield~-1+rep+gen, random=~rep:Row+rep:Col,randomMatrices=lM, data=df_ext)
obj5$logL
obj5$ED

dev5 = -2*obj5$logL + Constant
dev5

# sum over rows is zero, as B is zero matrix, so doesn't
# depend on parameters a
A1 <- kronecker(diag(nrep), kronecker(diag(nrows), rep(1,ncols)))
B1 <- t(A1) %*% lZ[[3]]
dim(B1)
range(B1)

# sum over cols is zero, as B is zero matrix, so doesn't
# depend on parameters a
A2 <- kronecker(diag(nrep), kronecker(rep(1,nrows), diag(ncols)))
B2 <- t(A2) %*% lZ[[3]]
dim(B2)
range(B2)

# sum over rows within replicates all zero:
sum(abs(B1 %*% coef(obj5)$ProwxPcol)<1.0e-12)
nrep*nrows

# sum over cols within replicates all zero:
sum(abs(B2 %*% coef(obj5)$ProwxPcol)<1.0e-12)
nrep*ncols

eff_Prowcol <- coef(obj5)$ProwxPcol
e <- (nrows-1)*(ncols-1)*c(1:nrep)
s <- c(1,e[-nrep]+1)
M <- NULL
for (r in 1:nrep)
{
  cat("rep ", r, ": \n")
  eff <- U12sc %*% eff_Prowcol[s[r]:e[r]]
  Mrep <- matrix(data = eff,nrow=nrows,ncol=ncols,byrow=TRUE)
  print(range(rowSums(Mrep)))
  print(range(colSums(Mrep)))
  M <- cbind(M, Mrep)
}

#pdf("Smooth_interaction_Gilmour_data.pdf",width=10,height=7)
# total numer of columns is 3*5 = 15
col = 1:(nrep*ncols)
row = 1:nrows
image.plot(col,row,t(M))
axis(side=3,at=c(3,8,13),labels=c('rep1','rep2','rep3'))
abline(v=5.5,lw=2.5)
abline(v=10.5,lw=2.5)
#dev.off()


