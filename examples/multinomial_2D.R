library(tidyverse)
library(LMMsolver)
library(ggplot2)
library(spam)

set.seed(1234)
n <- 250

nc <- 3
nc + 1    # total number of categories

mu1 <- c(0.3, 0.5, 0.7)
mu2 <- c(0.5, 0.4, 0.8)

sc <- c(4, 6, 4)

sim_non_linear <- function(x, mu1, mu2, sc) {
  z <- 2*exp(-4*((x[1]-mu1)^2+(x[2]-mu2)^2))-1
  return(sc*z)
}

sim_fun <- function(x, mu1, mu2, sc) {
  nCat <- length(mu1)
  z <- rep(NA, nCat)
  for (i in seq_len(nCat)) {
    z[i] <- sim_non_linear(x, mu1[i], mu2[i], sc[i])
  }
  fam <- multinomial()
  hz <- fam$linkinv(z)
  d <- c(hz, 1-sum(hz))
  return(d)
}

nobs <- 250
x1 <- runif(nobs)
x2 <- runif(nobs)
x <- data.frame(x1=x1,x2=x2)
prob <- t(apply(X=x, MARGIN=1, FUN=function(x) {sim_fun(x, mu1, mu2, sc)}))

sz <- 1 #sample(10:30, size=n, replace = TRUE)
M <- cbind(prob,sz)
multiNom <- t(apply(M, MARGIN=1, FUN=
                      function(x) {
                        rmultinom(n=1, size=x[nc+2], prob=x[1:(nc+1)])
                      } ))

colNames <- paste0(LETTERS[1:(nc+1)])
colnames(multiNom) <- colNames

dat <- data.frame(x1,x2,multiNom)
head(dat)

dat$category <-as.factor(LETTERS[apply(multiNom,MARGIN=1,FUN=function(x){which(x==1)})])
ggplot(dat,aes(x=x1,y=x2,color=category)) + geom_point()


x <- expand.grid(x1=seq(0,1,length=101),
                 x2=seq(0,1,length=101))

prob_true <- t(apply(X=x, MARGIN=1, FUN=function(x) {sim_fun(x, mu1, mu2, sc)}))

dat_true <- data.frame(x, prob_true)
colnames(dat_true) <- c("x1","x2",LETTERS[1:4])

dat_long <- gather(dat_true, category, prob, A:D, factor_key=TRUE)

ggplot(dat_long) +
  geom_tile(aes(x = x1, y = x2, fill = prob)) +
  scale_fill_viridis_c(option = "magma", limits=c(0,1)) +
  facet_wrap(vars(category), ncol = 2) + ggtitle("true")


obj <- LMMsolve(cbind(A,B,C,D)~1,spline=~spl2D(x1=x1,x2=x2, nseg=c(20,20),
                                               x1lim=c(0,1),x2lim=c(0,1)),
                data=dat, family=multinomial(),trace=TRUE)
summary(obj)
newdat <- expand.grid(x1=seq(0,1,length=101),
                 x2=seq(0,1,length=101))
pred <- predict(obj, newdata=newdat)

ggplot(pred) +
  geom_tile(aes(x = x1, y = x2, fill = ypred)) +
  scale_fill_viridis_c(option = "magma", limits=c(0,1)) +
  facet_wrap(vars(category), ncol = 2) + ggtitle("fitted")



