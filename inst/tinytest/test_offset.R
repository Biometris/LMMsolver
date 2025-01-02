set.seed(1234)

lambda <- 5000

nr <- c(0:3)
d <- 10^(-nr)
cntM <- sapply(d, FUN = function(x){ rpois(nrep,lambda*x)})
cnt <- as.vector(cntM)

dat <- data.frame(cnt = cnt, dilution = d)

# offset numeric
obj1 <- LMMsolve(cnt ~ 1, family = poisson(),
                offset = log(dat$dilution),
                data = dat)

mu <- coef(obj)$'(Intercept)'  # MB
lambda <- as.numeric(exp(mu))

# analytic solution
ML_lambda <- sum(dat$cnt)/sum(dat$dilution)

expect_equal(lambda, ML_lambda)

# offset is name of column in data frame
dat$log_dilution <- log(dat$dilution)
obj2 <- LMMsolve(cnt ~ 1, family = poisson(),
                offset = "log_dilution",
                data = dat)
mu <- coef(obj)$'(Intercept)'  #
lambda <- as.numeric(exp(mu))
expect_equal(lambda, ML_lambda)

expect_error(LMMsolve(cnt ~ 1, family = poisson(),
                 offset = "colname",
                 data = dat),
             "offset colname not defined in the data.")


