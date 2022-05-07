# Martin Boer, Biometris, WUR
#
# small example unbalanced data from
# Cunningham and Henderson (1968) Biometrics, reanalyzed using REML in
# Patterson and Thompson (1971) Biometrika:
#
#
rm(list=ls())
library(asreml)
library(LMMsolver)
library(mgcv)
library(nlme)
library(dplyr)

dat <- read.csv("cunningham_henderson_ex.csv", stringsAsFactors = TRUE)
head(dat)

# unbalanced data:
dat %>% group_by(treatment, block) %>% tally()

# asreml for comparison with LMMsolver
obj1 <- asreml(y~treatment, random = ~block, data=dat)
obj1$loglik
summary(obj1)$varcomp

coef(obj1, list=TRUE)

# LMMsolver:
obj2 <- LMMsolve(y~treatment, random=~block, data=dat)
obj2$logL
summary(obj2, which = "variances")

# effective dimension:
obj2$ED

coef(obj2)

# nlme4:
obj3 <- lme(y~treatment,data=dat,~1|block)

# mgcv:
obj4 <- gam(y ~ treatment + s(block,bs="re"), data=dat, method="REML")
summary(obj4)

#
# Compare the logL for the four models (asreml, LMMsolver, nlme, mgcv) used,
# should be all equal:
#
# direct comparison between asreml and LMMsolver, asreml doesn't include constant:
obj1$loglik
logLik(obj2, includeConstant = FALSE)

Constant_logL <- obj2$constantREML

round(obj1$loglik + Constant_logL,3)
round(obj2$logL + Constant_logL, 3)
round(obj3$logLik, 3)
round(as.numeric(-obj4$gcv.ubre),3)





