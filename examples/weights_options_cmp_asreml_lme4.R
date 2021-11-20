library(LMMsolver)
library(asreml)
library(lme4)

# set seed
set.seed(1234)

# number of treatments
nTrt <- 100

# residual variances and for treatment:
sigma2e <- 1.0
sigma2Trt <- 5.0

# block 1, two samples per treatment
nS1 <- rep(2, nTrt)

TrtEff <- rnorm(n=nTrt, sd=sqrt(sigma2Trt))
eps1 <- rnorm(n=nTrt, sd=sqrt(sigma2e/nS1))
y1 <- TrtEff + eps1

# block 2, 5 samples per treatments
nS2 <- rep(5, nTrt)

eps2 <- rnorm(n=nTrt, sd=sqrt(sigma2e/nS2))
y2 <- TrtEff + eps2

# combine the two blocks:
dat <- data.frame(block = rep(c("B1","B2"), each=nTrt),
                  trt = rep(paste0("Trt",1:nTrt), times=2),
                  y = c(y1,y2), stringsAsFactors = TRUE)

head(dat)

# some missing values for y:
dat$y[sample(x=c(1:nrow(dat)), size=4)] <- NA

# weight as number of samples, more samples, more weight (less penalty):
dat$w <- c(nS1, nS2)

obj1 <- asreml(fixed = y~block,
               random=~trt,
               weights = "w",
               data = dat)
summary(obj1)$varcomp

obj2 <- LMMsolve(fixed = y~block,
                random=~trt,
                weights = dat$w,
                data = dat,
                trace = FALSE,
                tolerance = 1.0e-10)
summary(obj2, which='variance')

obj3 <- lmer(y~block+(1|trt), weights = dat$w, data = dat)

# compare asreml-R with LMMsolver:
obj1$loglik
logLik(obj2, includeConstant = FALSE)

# compare LMMsolve with LME4:
logLik(obj2, includeConstant = TRUE)
logLik(obj3)


