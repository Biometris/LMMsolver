# Martin Boer, Biometris, WUR
#
# balanced one way ANOVA and the connection with effective dimensions
#
rm(list=ls())
library(asreml)
library(LMMsolver)
library(dplyr)

# a genotypes, n times replicated
a <- 10
n <- 3

set.seed(1234)

sigma2g <- 1.0
sigma2e <- 1.0

geno_eff  =rnorm(a,sd=sqrt(sigma2g))
geno_name = paste0("g",1:a)

dat <- data.frame(geno=rep(geno_name, each=n),
                  rep=rep(paste0("rep",1:n),times=a),
                  y = rep(geno_eff,each=n) + rnorm(a*n,sd=sigma2e),
                  stringsAsFactors = TRUE)

# asreml for comparison with LMMsolver
obj1 <- asreml(y~1, random = ~geno,data=dat, trace = FALSE)
obj1$loglik

# LMMsolver:
obj2 <- LMMsolve(y~1, random=~geno, data=dat, eps=1.0e-12)
obj2$logL

# effective dimension:
obj2$ED

# calculation of MSA and MSE by foot, see
# p. 70-71 Searle et al. Variance Components.
y <- dat$y
I_a <- diag(a)
I_n <- diag(n)
J_a <- (1/a) * matrix(data=1,ncol=a,nrow=a)
J_n <- (1/n) * matrix(data=1,ncol=n,nrow=n)
SSA <- as.numeric(t(y) %*% (I_a %x% J_n - J_a %x% J_n) %*% y)
SSE <- as.numeric(t(y) %*% (I_a %x% I_n - I_a %x% J_n) %*% y)
MSA <- SSA/(a-1)
MSE <- SSE/(a*(n-1))

sigma2e_est <- MSE
sigma2g_est <- (MSA-MSE)/n

obj2$sigma2e
sigma2e_est

obj2$tau2e
sigma2g_est

# effective dimensions:
Fval <- MSA/MSE  # see page 64,65 Searle et al. Variance components:
EDg <-  (a-1)   - (1/Fval) * (a-1)
EDr <- a*(n-1)  + (1/Fval) * (a-1)

obj2$ED
EDg
EDr

# calculate F-distribution and compare with anova analysis,
# see page 64 Searle et al.
pf(Fval, df1=a-1, df2=a*(n-1), lower.tail=FALSE)

lm.obj <- lm(y~geno, data=dat)
anova(lm.obj)

# for balanced one way ANOVA we can use effective dimension from mixed model
# analysis to calculate F statistic:
EDgMax <- a-1
EDrMin <- (a*n-1) - (a-1)
Fval2 <- (EDgMax)/(EDgMax - obj2$ED['geno'])
pf(Fval2, df1=EDgMax, df2=EDrMin, lower.tail=FALSE)



