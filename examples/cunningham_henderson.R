# Martin Boer, Biometris, WUR
#
# small example unbalanced data from
# Cunningham and Henderson (1968) Biometrics, reanalyzed using REML in
# Patterson and Thompson (1971) Biometrika:
#
rm(list=ls())
library(asreml)
library(LMMsolver)
library(dplyr)

dat <- read.csv("cunningham_henderson_ex.csv")
head(dat)

# unbalanced data:
dat %>% group_by(treatment, block) %>% tally()

# asreml for comparison with LMMsolver
obj1 <- asreml(y~treatment, random = ~block,data=dat)
obj1$loglik

coef(obj1, list=TRUE)

# LMMsolver:
obj2 <- LMMsolve(y~-1+treatment, random=~block, data=dat)
obj2$logL

# effective dimension:
obj2$ED

coef(obj2)

