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

D1 <- diff(diag(nrows), diff=1)
DtD1 <- crossprod(D1)
U1 <- eigen(DtD1)$vectors[,1:(nrows-1)]
d1 <- eigen(DtD1)$values[1:(nrows-1)]

U1sc <- U1 %*% diag(1/sqrt(d1))
tmp1 <- kronecker(U1sc, rep(1,ncols))
Z3 = kronecker(diag(nrep),tmp1)

df_ext = cbind(df,Z3)

lM <- list(Prow=c(8:70))
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

D2 <- diff(diag(ncols), diff=1)
DtD2 <- crossprod(D2)
U2 <- eigen(DtD2)$vectors[,1:(ncols-1)]
d2 <- eigen(DtD2)$values[1:(ncols-1)]

U2sc <- U2 %*% diag(1/sqrt(d2))
tmp2 <- kronecker(rep(1,nrows), U2sc)
Z4 = kronecker(diag(nrep),tmp2)

df_ext <- cbind(df_ext, Z4)
dim(df_ext)
lM <- list(Prow=c(8:70),Pcol=c(71:82))
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

U12sc = kronecker(U1 %*% diag(1/sqrt(d1)), U2 %*% diag(1/sqrt(d2)))

Z5 = kronecker(diag(nrep), U12sc)

dim(Z5)
df_ext <- cbind(df_ext, Z5)
dim(df_ext)

lM <- list(Prow=c(8:70),Pcol=c(71:82),ProwxPcol=c(83:334))
obj5 = LMMsolve(fixed=yield~-1+rep+gen, random=~rep:Row+rep:Col,randomMatrices=lM, data=df_ext)
obj5$logL
obj5$ED

dev5 = -2*obj5$logL + Constant
dev5

# sum over rows is zero, as B is zero matrix, so doesn't
# depend on parameters a
A1 <- kronecker(diag(nrep), kronecker(diag(nrows), rep(1,ncols)))
B1 <- t(A1) %*% Z5
dim(B1)
range(B1)

# sum over cols is zero, as B is zero matrix, so doesn't
# depend on parameters a
A2 <- kronecker(diag(nrep), kronecker(rep(1,nrows), diag(ncols)))
B2 <- t(A2) %*% Z5
dim(B2)
range(B2)

sum(abs(B1 %*% coef(obj5)$ProwxPcol)<1.0e-12)
sum(abs(B2 %*% coef(obj5)$ProwxPcol)<1.0e-12)

eff_Prowcol <- coef(obj5)$ProwxPcol
rep1 <- U12sc %*% eff_Prowcol[1:84]
rep2 <- U12sc %*% eff_Prowcol[85:168]
rep3 <- U12sc %*% eff_Prowcol[169:252]

Mrep1 <- matrix(data=rep1,nrow=22,ncol=5,byrow=TRUE)
rowSums(Mrep1)
colSums(Mrep1)

Mrep2 <- matrix(data=rep2,nrow=22,ncol=5,byrow=TRUE)
rowSums(Mrep2)
colSums(Mrep2)

Mrep3 <- matrix(data=rep3,nrow=22,ncol=5,byrow=TRUE)
rowSums(Mrep3)
colSums(Mrep3)

M <- cbind(Mrep1, Mrep2, Mrep3)

#pdf("Smooth_interaction_Gilmour_data.pdf",width=10,height=7)
col = 1:15
row = 1:22
image.plot(col,row,t(M))
axis(side=3,at=c(3,8,13),labels=c('rep1','rep2','rep3'))
abline(v=5.5,lw=2.5)
abline(v=10.5,lw=2.5)
#dev.off()

# sum equal to zero, follows by definition from U2sc
eff_Prow <- coef(obj5)$Prow
eff_Pcol <- coef(obj5)$Pcol

col_eff_rep1 = U2sc %*% eff_Pcol[1:4]
col_eff_rep2 = U2sc %*% eff_Pcol[5:8]
col_eff_rep3 = U2sc %*% eff_Pcol[9:12]

row_eff_rep1 <- U1sc %*% eff_Prow[1:21]
row_eff_rep2 <- U1sc %*% eff_Prow[22:42]
row_eff_rep3 <- U1sc %*% eff_Prow[43:63]

row_eff_rep1_wh = U1sc %*% coef(obj4)$Prow[1:21]
col_eff_rep1_wh = U2sc %*% coef(obj4)$Pcol[1:4]

plot(x=1:22, y=row_eff_rep1,type='b',col='red',ylim=c(-60,60),lwd=2.5)
points(x=1:22, y=row_eff_rep1_wh,type='b',col='green',lwd=2.5)
for (i in 1:5) {
  points(x=1:22,y=row_eff_rep1 + Mrep1[,i],type='b')
}

plot(x=1:5, y=col_eff_rep1,type='b',col='red',ylim=c(-150,150))
points(x=1:5, y=col_eff_rep1_wh,type='b',col='green')

for (i in 1:22) {
  points(x=1:5,y=col_eff_rep1 + Mrep1[i,])
}

