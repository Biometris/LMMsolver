#
# Martin Boer, Lettuce example with population mean zero,
#              plus some simulated data.
#
rm(list = ls())
library(fields)
library(spam)
library(LMMsolver)
library(dplyr)
library(ggplot2)

dat = read.csv('lettuce_plus_harvest.csv',header=TRUE, stringsAsFactors = FALSE)

head(dat)
dat$x = dat$harvestday/365
labels <- unique(dat$geno)
dat$geno <- factor(dat$geno,levels=labels)
dat$env <- as.factor(dat$env)
dat$g_nr <- as.numeric(dat$geno)

Ngeno <- nlevels(dat$geno)
Nenv <- nlevels(dat$env)

# just some extra simulation to add slope and some missing values..
set.seed(1234)
slope_G <- rnorm(Ngeno,sd=0.25)
slope_G_ext <- rep(slope_G,each=Nenv)
dat$ysim = dat$nitrate + dat$x * slope_G_ext

# one value missing...
ndx <- sample(x=nrow(dat),size=3)
dat <- dat[-ndx,]

obj0 <- LMMsolve(ysim~env,random=~geno, data=dat)
obj0$ED

# a very simple simple B-spline basis, helps that
# data can be in any order
knots2 = PsplinesKnots(1, Ngeno, 1, Ngeno-1)
B2 <- Bsplines(knots2, dat$g_nr)

# define matrix orthogonal to constant
J = matrix(data=1,nrow=Ngeno, ncol=Ngeno)
K = diag(Ngeno) - (1/Ngeno) * J
U2sc = eigen(K)$vectors[,-Ngeno]
t(U2sc) %*% rep(1,Ngeno)

lZ <- list()
lZ[[1]] = B2 %*% U2sc
Z <- do.call("cbind", lZ)
dat_ext = cbind(dat, Z)

lM <- ndxMatrix(dat, lZ, c("geno"))
obj1 <- LMMsolve(ysim~env, randomMatrices=lM, data=dat_ext)

obj0$ED
obj1$ED

obj0$logL
obj1$logL

degr1 = 3
pord1 = 2
xmin1 <- 0.2
xmax1 <- 1.5
nseg1 <- 10

knots1 = PsplinesKnots(xmin1, xmax1, degr1, nseg1)
B1  <- Bsplines(knots1, dat$x)
q1 <- ncol(B1)

U1sc <- calcUsc(q1, pord1)

lZ <- list()
lZ[[1]] = B1 %*% U1sc  #env
lZ[[2]] = B2 %*% U2sc  #geno
lZ[[3]] = dat$x*(B2 %*% U2sc) # x.geno
lZ[[4]] = RowKronecker(B1,B2) %*% kronecker(U1sc,U2sc)

Z <- do.call("cbind", lZ)
dat_ext = cbind(dat, Z)

lM <- ndxMatrix(dat, lZ, c("f(z)","g","G.x","f_g(z)"))
obj5 <- LMMsolve(ysim~x, randomMatrices=lM, data=dat_ext)
round(obj5$ED, 2)

x0 <- seq(xmin1,xmax1,by=0.01)
B1grid <- Bsplines(knots1,x0)
mu <- coef(obj5)$'(Intercept)'
beta <- coef(obj5)$x
ypredEnv <- mu + beta*x0 + (B1grid %*% U1sc) %*% coef(obj5)$'f(z)'
ypredE.df <- data.frame(x=x0, y=ypredEnv)
p <- ggplot(dat, aes(x=x, y=nitrate, color=geno)) + geom_point()
p + geom_line(ypredE.df, mapping=aes(x=x,y=y),col='black',size=1.5)

G_eff <- as.vector(U2sc %*% coef(obj5)$G)
G.x <- as.vector(U2sc %*% coef(obj5)$G.x)
ypredG <- as.vector(kronecker(rep(1,length(x0)), G_eff))
ypredGx <- as.vector(kronecker(x0, G.x))
ypredGxE <- ypredG + ypredGx +
            kronecker(B1grid,diag(Ngeno)) %*% kronecker(U1sc,U2sc) %*% coef(obj5)$'f_g(z)'

M <- matrix(data=ypredGxE,nrow=length(x0),ncol=Ngeno,byrow=TRUE)
range(rowSums(M))

x1_ext <- kronecker(x0, rep(1,Ngeno))
x2_ext <- kronecker(rep(1,length(x0)), c(1:Ngeno))
pred <- data.frame(x1=x1_ext,x2=x2_ext, ypredGxE)

ggplot(pred)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=ypredGxE))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
  ggtitle("Genotypic deviations from mean") + xlab("z") + ylab("genotype") +
  theme(plot.title = element_text(hjust = 0.5))

