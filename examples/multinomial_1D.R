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
head(dat)

sRows <- rowSums(multiNom)
fr <- multiNom/sRows
dat_fr <- data.frame(x, fr)

obj <- LMMsolve(fixed = cbind(A,B,C,D) ~ 1,
                spline = ~spl1D(x, nseg = 17, xlim=c(0,1), scaleX=FALSE),
                data=dat, family = multinomial())
summary(obj)
coef(obj)
deviance(obj)

# make predictions
x0 <- seq(0, 1, by = 0.01)
newdat <- data.frame(x=x0)
pred <- predict(obj, newdata=newdat)
colnames(pred) <- c("x", "category", "y")
head(pred)

prob_true <- t(sapply(X=x0, FUN=function(x) { sim_fun(x, mu, sc)}))
colnames(prob_true) <- colNames
nCatTot <- length(obj$respVar)
df_true <- data.frame(x0, prob_true)
prob_true_lf <- data.frame(category = rep(colNames, each=length(x0)),
                           x=rep(x0,times=nCatTot),
                           y = as.vector(prob_true))

dat_fr_lf <- data.frame(category = rep(colNames, each=nrow(fr)),
                     x=rep(x,times=nCatTot),
                     y = as.vector(fr))

p1 <- ggplot(prob_true_lf, aes(x = x, y=y,color=category)) +
  geom_line(linetype='dashed') +
  geom_line(data=pred) +
  geom_point(data=dat_fr_lf)
p1

