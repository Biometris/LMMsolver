library(LMMsolver)
library(ggplot2)
library(spam)

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
  z <- rep(NA, 3)
  for (i in 1:3) {
    z[i] <- sim_non_linear(x, mu[i], sc[i])
  }
  hz <- LMMsolver:::glogit(z)
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
colnames(multiNom) <- paste0(LETTERS[1:(nc+1)])

dat <- data.frame(x, multiNom)
head(dat)

sRows <- rowSums(multiNom)
fr <- multiNom/sRows
dat_fr <- data.frame(x, fr)

obj <- LMMsolve(fixed = cbind(A,B,C,D) ~ 1,
                spline = ~spl1D(x, nseg = 17, xlim=c(0,1), scaleX=FALSE),
                data=dat, family = multinomial())
# seems consistent now:
summary(obj)

# not correct reference to columns in MME:
obj$ndxCoefficients

x0 <- seq(0, 1, by=0.01)
X0 <- cbind(1, x0) %x% diag(nc)
knots <- obj$splRes[[1]]$knots[[1]]
B0 <- LMMsolver:::Bsplines(knots, x0)
Z0 <- B0 %x% diag(nc)
beta <- obj$coefMME
eta0 <- cbind(X0, Z0) %*% beta
etaM <- matrix(data=eta0,nrow = length(x0), ncol= nc, byrow=TRUE)
pi_est <- t(apply(etaM, MARGIN=1, FUN=LMMsolver:::glogit))
pi_est <- cbind(pi_est, 1.0 - rowSums(pi_est))
colnames(pi_est) <- obj$respVar
pred <- data.frame(x=x0, pi_est)

prob_true <- t(sapply(X=x0, FUN=function(x) { sim_fun(x, mu, sc)}))
colnames(prob_true) <- obj$respVar
df_true <- data.frame(x0, prob_true)

p1 <- ggplot(df_true, aes(x = x0)) +
  geom_line(aes(y = A), color = "darkred", linetype='dashed') +
  geom_line(aes(y = B), color = "steelblue", linetype='dashed') +
  geom_line(aes(y = C), color = "green", linetype='dashed') +
  geom_line(aes(y = D), color = "black", linetype='dashed') +
  geom_line(data=pred, aes(x=x, y = A), color = "darkred") +
  geom_line(data=pred, aes(y = B), color = "steelblue") +
  geom_line(data=pred, aes(y = C), color = "green") +
  geom_line(data=pred, aes(y = D), color = "black") +
  geom_point(data=dat_fr, aes(x=x, y= A), color = "darkred") +
  geom_point(data=dat_fr, aes(x=x, y= B), color = "steelblue") +
  geom_point(data=dat_fr, aes(x=x, y= C), color = "green") +
  geom_point(data=dat_fr, aes(x=x, y= D), color = "black") +
  xlab("x") + ylab("prob") +
  geom_hline(yintercept=1.0, linetype='dashed')
p1



