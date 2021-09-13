# Martin Boer, Biometris

# using asreml4 for comparison with LMMsolver:
library(LMMsolver)
library(asreml)
library(dplyr)

df <- read.csv("multipopQTL.csv", stringsAsFactors = TRUE)
head(df)
str(df)

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
obj0 <- LMMsolve(fixed=pheno~cross, residualterm='cross', data=df, eps=1.0e-8,
                          monitor=TRUE, display= FALSE)
# check with asreml:
obj0$logL
obj0.asr$loglik

# effective degrees of freedom
df %>% group_by(cross) %>% tally()
obj0$ED

# include QTL, using LMMsolve
lM <- list(QTL=c(3:5))
obj1 <- LMMsolve(fixed=pheno~cross, group=lM,residualterm='cross', data=df, eps=1.0e-8,
                 monitor=TRUE, display= FALSE)

# check with asreml:
obj1$logL
obj1.asr$loglik

# effective degrees of freedom: maximum effective dimension for QTL is 2,
# because of sum to one constraint:
obj1$ED

# coefficients for QTL-effect:
coefficients(obj1.asr,list=TRUE)
coef(obj1)

coefficients(obj1.asr,list=TRUE)
coef(obj1)

# sum of effects equal to zero:
sum(coef(obj1)$QTL)

