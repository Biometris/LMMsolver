# Martin Boer, Biometris

library(LMMsolver)
library(asreml)
library(dplyr)

data(multipop)

df <- multipop

# NULL MODEL, no marker:
obj0.asr = asreml(fixed = pheno~cross,
                  residual = ~dsum(~units|cross),
                  data=df, trace=TRUE)
summary(obj0.asr)$varcomp

# alternative model, using 1 marker...
Lgrp <- list(QTL=c(3:5))
obj1.asr = asreml(fixed = pheno~cross,random=~grp(QTL),
                  residual = ~dsum(~units|cross),
                  group = Lgrp[],
                  data=df,trace=TRUE)
summary(obj1.asr)$varcomp
dev = 2.0*obj1.asr$log - 2.0*obj0.asr$log
dev
minlog10p = -log10(0.5*pchisq(dev,1,lower.tail=FALSE))
minlog10p

# NULL mode, using LMMsolver:
obj0 <- LMMsolve(fixed=pheno~cross, residual=~cross, data=df, tolerance=1.0e-8,
                          trace=TRUE)
# check with asreml:
obj0$logL
obj0.asr$loglik
obj0$theta

# include QTL, using LMMsolve
lM <- list(QTL=c(3:5))
obj1 <- LMMsolve(fixed=pheno~cross, group=lM,
                 random = ~grp(QTL),
                 residual=~cross,
                 data=df, tolerance=1.0e-6,
                 trace=TRUE)

# check with asreml:
obj1$logL
obj1.asr$loglik

# strong QTL effect, eff dimension close to two:
summary(obj1)

# coefficients for QTL-effect asreml:
coefficients(obj1.asr, list=TRUE)

# coefficients with LMMsolver
coef(obj1)

# sum of effects equal to zero:
sum(coef(obj1)$QTL)

obj0$ED

obj1$ED
obj1$EDmax
obj1$EDnominal

summary(obj1.asr, coef=TRUE)

theta <- obj1$theta
listC <- list()
X = model.matrix(~cross, data=df)
Z = as.spam(as.matrix(df[,3:5]))
U = as.spam(cbind(X,Z))
Rinv <- LMMsolver:::constructRinv(df, residual=~cross,weights=1.0)
Ut = t(U)
Ginv = diag.spam(c(0,0,1,1,1))
listC <-list()
listC[[1]] = Ginv
listC[[2]] = Ut %*% Rinv[[1]] %*% U
listC[[3]] = Ut %*% Rinv[[2]] %*% U
objAD <- LMMsolver:::ADchol(listC)

ED <- theta * LMMsolver:::dlogdet(objAD, theta)
c(3, 100, 80) - ED

C <- LMMsolver:::linearSum(theta,listC)
cholC <- chol(C)
A <- LMMsolver:::DerivCholesky(cholC,objAD)
A
sqrt(diag(A))

pred1 <- predict(obj1.asr, classify='grp(QTL)',sed=TRUE)
