set.seed(1234)

set.seed(1234)
n <- 250

nc <- 3
nc + 1    # total number of categories

mu <- c(0.1, 0.3, 0.5)
sc <- c(4, 5, 4)

sim_non_linear <- function(x, mu, sc) {
  z <- 2*exp(-20*(x-mu)^2)-1
  return(sc*z)
}

sim_fun <- function(x, mu, sc) {
  nCat <- length(mu)
  z <- rep(NA, nCat)
  for (i in seq_len(nCat)) {
    z[i] <- sim_non_linear(x, mu[i], sc[i])
  }
  fam <- multinomial()
  hz <- fam$linkinv(z)
  d <- c(hz, 1-sum(hz))
  return(d)
}

x <- runif(n, 0, 1)   #seq(0, 1, length=n)
prob <- t(sapply(X=x, FUN=function(x) {sim_fun(x, mu, sc)}))
range(rowSums(prob))

sz <- sample(10:30, size=n, replace = TRUE)
M <- cbind(prob,sz)
multiNom <- t(apply(M, MARGIN=1, FUN=
                      function(x) {
                        rmultinom(n=1, size=x[nc+2], prob=x[1:(nc+1)])
                      } ))

colNames <- paste0(LETTERS[1:(nc+1)])
colnames(multiNom) <- colNames

dat <- data.frame(x, multiNom)

mod1 <- LMMsolve(fixed = cbind(A,B,C,D) ~ 1,
                spline = ~spl1D(x, nseg = 17, xlim=c(0,1), scaleX=FALSE),
                data=dat, family = multinomial())
expect_equivalent_to_reference(mod1, "multinomial1")

dat2 <- dat
dat2$A[1] <- NA

expect_error(LMMsolve(fixed = cbind(A,B,C,D) ~ 1,
                 spline = ~spl1D(x, nseg = 17, xlim=c(0,1), scaleX=FALSE),
                 data=dat2, family = multinomial()),
             "family multinomial : NA's in response variable.")

dat3 <- dat
dat3$A[1] <- -1

expect_error(LMMsolve(fixed = cbind(A,B,C,D) ~ 1,
                      spline = ~spl1D(x, nseg = 17, xlim=c(0,1), scaleX=FALSE),
                      data=dat3, family = multinomial()),
             "family multinomial : negative values in response variable.")

dat4 <- dat
dat4$A[1] <- "A"

expect_error(LMMsolve(fixed = cbind(A,B,C,D) ~ 1,
                      spline = ~spl1D(x, nseg = 17, xlim=c(0,1), scaleX=FALSE),
                      data=dat4, family = multinomial()),
             "response cbind(A, B, C, D) should be numeric.", fixed=TRUE)

expect_error(LMMsolve(fixed = cbind(A,B) ~ 1,
                      spline = ~spl1D(x, nseg = 17, xlim=c(0,1), scaleX=FALSE),
                      data=dat, family = multinomial()),
             "family multinomial two categories, use binomial family.", fixed=TRUE)





