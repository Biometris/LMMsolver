# Martin Boer, Biometris
library(asreml)
library(LMMsolver)

# two crosses, AxB an AxC
df <- read.csv("Example3QTLs.csv",stringsAsFactors = TRUE)
str(df)

Lgrp <- list(QTL1=c(3:5), QTL2=c(6:8), QTL3=c(9:11))

obj1.asr <- asreml(fixed = pheno~cross,random=~grp(QTL1)+grp(QTL2)+grp(QTL3),
                  residual = ~dsum(~units|cross),
                  group = Lgrp,
                  data=df,trace=TRUE)

obj1.LMM <- LMMsolve(fixed=pheno~cross, random=~grp(QTL1)+grp(QTL2)+grp(QTL3),
                 group = Lgrp,
                 residual=~cross, data=df,trace=TRUE)
summary(obj1.LMM)

coef(obj1.LMM)

# only small difference, ok
obj1.asr$loglik
obj1.LMM$logL
obj1.LMM$logL - obj1.asr$loglik

