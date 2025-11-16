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

expect_error(obtainSmoothTrend(mod1,grid=100),
             "For multinomial response use predict function.")

newdat <- data.frame(x=seq(0,1,by=0.01))
expect_error(predict(mod1, newdata=newdat, se.fit=TRUE),
             "se.fit=TRUE not implemented yet for multinomial.")

# example with factor
dat$Factor <- as.factor(ifelse(dat$x < 0.5, "fA","fB"))
mod2 <- LMMsolve(fixed = cbind(A,B,C,D) ~ Factor,
                 spline = ~spl1D(x, nseg = 17, xlim=c(0,1), scaleX=FALSE),
                 data=dat, family = multinomial())
newdat$x <- as.factor(ifelse(newdat$x < 0.5, "fA","fB"))
expect_error(predict(mod2, newdata=newdat),
             "use of factors not implemented yet for multinomial.")

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


dat <- data.frame(y=c(1,2,3,2,1,2))

expect_error(LMMsolve(y~1, family = multinomial(), data=dat),
             "family multinomial : response should be a matrix or a factor.", fixed=TRUE)



set.seed(1234)
sz <- 10
n <- 100
multiNom <- t(rmultinom(n=n,size=sz,prob=c(0.1,0.2,0.3,0.4)))
colnames(multiNom) <- LETTERS[1:4]
dat <- data.frame(multiNom)

obj <- LMMsolve(cbind(A,B,C,D)~1, data=dat, family=multinomial())
cf <- coef(obj)$'(Intercept)'
fam <- multinomial()
pr <- fam$linkinv(cf)
est_D <- sum(1-sum(pr))
prob <- as.numeric(c(pr, est_D))
prob_ML <- as.numeric(colSums(multiNom)/(sz*n))
# there are small differences, LMMsolver uses a small probability
# that the catagories have wrong labels, to keep algorithm stable
expect_equal(prob, prob_ML, tolerance = 1.0e-6)

# check for weight argument
expect_error(LMMsolve(fixed = cbind(A,B,C,D) ~ 1,
         spline = ~spl1D(x, nseg = 17, xlim=c(0,1), scaleX=FALSE),
         data=dat, family = multinomial(), weights=c(1,nrow(dat))),
         "family multinomial : weights option cannot be used.")

