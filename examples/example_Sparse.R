library(LMMsolver)
library(nlme)

set.seed(1234)

# number of levels
K <- 1000
r <- 5
n <- K*r

s <- rnorm(K, sd=1.0)
e <- rnorm(n, sd=5.0)

y = rep(s,r) + e

df <- data.frame(y,block= as.factor(rep(paste0("blk",1:K), r)),
                   rep  = as.factor(rep(paste0("rep",1:r), each=K)))
head(df)
str(df)

# This model doesn't work for K > 1000, problem with intercept:
#obj0 <- LMMsolve(y~1,random=~block, data=df, trace=TRUE)

# This model is ok, with replicate fixed:
obj1 <- LMMsolve(y~1,random=~block, data=df, trace=TRUE)
obj1$ED
# nlme4:
obj2 <- lme(y~1,data=df,~1|block)

obj1$logL

p <- 1
N <- nrow(df)

Constant_logL <- -0.5*log(2*pi)*(N-p)

logL1 <- obj1$logL + Constant_logL
logL2 <- obj2$logLik
logL1 - logL2


