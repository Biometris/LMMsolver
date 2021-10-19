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

