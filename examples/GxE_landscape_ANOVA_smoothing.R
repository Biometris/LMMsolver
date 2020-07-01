#' ---
#' title: GxE analysing using ANOVA smoothing
#' author: Martin Boer
#' ---

rm(list=ls())

library(ggplot2)
library(sampling)
library(splines)
library(LMMsolver)

set.seed(1234)

# residual error
sigma2 <- 4.0

# simulated surface a rotated ellipse
f <- function(x1,x2, theta) {
  a = 0.2
  b = 0.3
  x1c = 0.4
  x2c = 0.6
  x1r <- cos(theta)*(x1-x1c) - sin(theta)*(x2-x2c)
  x2r <- sin(theta)*(x1-x1c) + cos(theta)*(x2-x2c)
  20 - (x1r^2/a^2 + x2r^2/b^2)
}

# rotation of the simulated ellipse surface:
degrees_rotation = 30

x1grid <- (c(1:100)-0.5)/100
x2grid <- (c(1:100)-0.5)/100
grid <- expand.grid(x1grid, x2grid)
names(grid) <- c("x1","x2")

# 12 environments:
# 25 genotypes:
nEnv = 12
nGeno = 25

x1f <- findInterval(x1grid,seq(from=0,to=1,length=nEnv+1))
ndx1 <- strata(data.frame(x1f),stratanames="x1f",size=rep(1,times=nEnv),method='srswor')$ID_unit

x2f <- findInterval(x1grid,seq(from=0,to=1,length=nGeno+1))
ndx2 <- strata(data.frame(x2f),stratanames="x2f",size=rep(1,times=nGeno),method='srswor')$ID_unit

grid$ysim <- f(grid$x1, grid$x2, theta=(degrees_rotation/180)*pi)

#env_cov <- runif(nEnv,min=0.02,max=0.98)
#geno_cov <- runif(nGeno,min=0.02,max=0.98)
#env_cov <- seq(from=0.05, to=0.95, length=nEnv)
#geno_cov <- seq(from=0.05, to=0.95, length=nGeno)

env_cov <- x1grid[ndx1]
geno_cov <- x2grid[ndx2]

env = paste0("env",formatC(1:nEnv,width=2,flag="0"))
gen = paste0("gen",formatC(1:nGeno,width=2,flag="0"))
env_df <- data.frame(env,z=env_cov)
geno_df <- data.frame(gen,x=geno_cov)

#
sim.df <- expand.grid(env_cov, geno_cov)
names(sim.df) <- c("x1","x2")
Nobs <- nrow(sim.df)
sim.df$ysim <- f(sim.df$x1, sim.df$x2, theta=(degrees_rotation/180)*pi)
e <- rnorm(Nobs, mean = 0, sd = sqrt(sigma2))
sim.df$y <- sim.df$ysim + e

ggplot(grid)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=ysim))+
  scale_fill_gradientn(name="ysim",colours= topo.colors(100)) +
#  scale_fill_gradient(low = "blue", high = "red",name="z") +
  #geom_vline(xintercept=s1bnd) +
  #geom_hline(yintercept=s1bnd) +
  geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("simulated") + xlab("z") + ylab("x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()

# row wise kronecker product
Rten2 <-
  function(X1,X2) {
    one.1 <- matrix(1,1,ncol(X1))
    one.2 <- matrix(1,1,ncol(X2))
    kronecker(X1,one.2)*kronecker(one.1,X2)
  }

degree = 3
k = degree+1
x1min = 0
x1max = 1
x2min = 0
x2max = 1
nseg1 = 10
nseg2 = 8
#nseg1 = 25
#nseg2 = 30

dx1 = (x1max - x1min) / nseg1
dx2 = (x2max - x2min) / nseg2

x1 = sim.df$x1
x2 = sim.df$x2

knots1 = seq(x1min - degree * dx1, x1max + degree * dx1, by = dx1)
knots2 = seq(x2min - degree * dx2, x2max + degree * dx2, by = dx2)

#knots1 <- c(rep(x1min,3),seq(x1min,x1max,by=dx1),rep(x1max,3))
#knots2 <- c(rep(x2min,3),seq(x2min,x2max,by=dx2),rep(x2max,3))

B1 = splineDesign(knots1, x1, derivs=rep(0,length(x1)), ord = degree+1)
q1 = ncol(B1)
q1

B2 = splineDesign(knots2, x2, derivs=rep(0,length(x2)), ord = degree+1)
q2 = ncol(B2)
q2

D1 <- diff(diag(q1), diff=1)
DtD1 <- crossprod(D1)
U1 <- eigen(DtD1)$vectors[,-q1]
d1 <- eigen(DtD1)$values[-q1]
U1sc <- U1 %*% diag(1/sqrt(d1))

D2 <- diff(diag(q2), diff=1)
DtD2 <- crossprod(D2)
U2 <- eigen(DtD2)$vectors[,-q2]
d2 <- eigen(DtD2)$values[-q2]
U2sc <- U2 %*% diag(1/sqrt(d2))

# J1 = matrix(data=1,ncol=q1,nrow=q1)
# K1 = diag(q1)- J1
# U1 <- eigen(K1)$vectors[,-q1]
# d1 <- eigen(K1)$values[-q1]
# U1sc <- U1 %*% diag(1/sqrt(d1))
#
# J2 = matrix(data=1,ncol=q2,nrow=q2)
# K2 = diag(q2)- J2
# U2 <- eigen(K2)$vectors[,-q2]
# d2 <- eigen(K2)$values[-q2]
# U2sc <- U2 %*% diag(1/sqrt(d2))

U12sc = kronecker(U1sc, U2sc)

df <- data.frame(gen=rep(gen,each=nEnv),
                 env=rep(env,times=nGeno),
                 y=sim.df$y, x1=sim.df$x1,x2=sim.df$x2,InterceptMB=1)
df_csv <- data.frame(gen=rep(gen,each=nEnv),
                env=rep(env,times=nGeno),
                y=sim.df$y, x=sim.df$x2,z=sim.df$x1)
head(df_csv)
write.csv(df_csv,"GxE_simulated_surface.csv",quote=FALSE,row.names=FALSE)

head(df)

lZ <- list()
lZ[[1]] = B1 %*% U1sc
lZ[[2]] = B2 %*% U2sc
Z <- do.call("cbind",lZ)

df_ext = cbind(df,Z)

dim <- sapply(lZ,ncol)
e <- cumsum(dim) + ncol(df)
s <- e - dim + 1

# model without interactions:
lM <- list(E=c(s[1]:e[1]), G=c(s[2]:e[2]))
obj0 = LMMsolve(fixed=y~-1+InterceptMB, random=NULL,randomMatrices=lM, data=df_ext)
obj0$logL
obj0$ED

lZ[[3]] = Rten2(B1,B2) %*% U12sc
Z <- do.call("cbind",lZ)
df_ext = cbind(df,Z)
dim <- sapply(lZ,ncol)
e <- cumsum(dim) + ncol(df)
s <- e - dim + 1

lM <- list(E=c(s[1]:e[1]), G=c(s[2]:e[2]), GxE=c(s[3]:e[3]))

obj1 = LMMsolve(fixed=y~-1+InterceptMB, random=NULL,randomMatrices=lM, data=df_ext)
obj1$logL
obj1$ED

theta1 <- U1sc %*% coef(obj1)$E
theta2 <- U2sc %*% coef(obj1)$G
sum(theta1)
sum(theta2)
theta3 <- U12sc %*% coef(obj1)$GxE

# marginal equal to zero:
M = matrix(data=theta3,nrow=q1,ncol=q2,byrow=TRUE)
rowSums(M)
colSums(M)

B1x = splineDesign(knots1, x1grid, derivs=rep(0,length(x1grid)), ord = degree+1)
B2x = splineDesign(knots2, x2grid, derivs=rep(0,length(x2grid)), ord = degree+1)

nGridext = 500
x1grid_ext = seq(x1min - degree * dx1, x1max + degree * dx1,length=nGridext)
x2grid_ext = seq(x2min - degree * dx2, x2max + degree * dx2,length=nGridext)
#x1grid_ext = seq(x1min, x1max,length=nGridext)
#x2grid_ext = seq(x2min, x2max,length=nGridext)

B1x_ext = splineDesign(knots1, x1grid_ext, derivs=rep(0,length(x1grid_ext)), ord = degree+1,outer.ok=TRUE)
B2x_ext = splineDesign(knots2, x2grid_ext, derivs=rep(0,length(x2grid_ext)), ord = degree+1,outer.ok=TRUE)

# marginals
sum(B1x_ext %*% theta1)
sum(B2x_ext %*% theta2)
plot(x=x1grid_ext,y=B1x_ext %*% theta1,type='l')
abline(h=0,col='red')
abline(v=x1min,lt=2)
abline(v=x1max,lt=2)

plot(x=x2grid_ext,y=B2x_ext %*% theta2,type='l')
abline(h=0,col='red')
abline(v=x1min,lt=2)
abline(v=x1max,lt=2)


ypred_ext <- kronecker(B1x_ext,B2x_ext) %*% theta3
M = matrix(ypred_ext,ncol=nGridext,nrow=nGridext)
range(rowSums(M))
range(colSums(M))

k1 <- length(x1grid_ext)
k2 <- length(x2grid_ext)
x1_ext <- kronecker(x1grid_ext, rep(1,k2))
x2_ext <- kronecker(rep(1,k1), x2grid_ext)

pred_ext <- data.frame(x1=x1_ext,x2=x2_ext,y=kronecker(B1x_ext,B2x_ext)%*%theta3)

ggplot(pred_ext)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=y))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
  #  scale_fill_gradient(low = "blue", high = "red",name="Fitted") +
  geom_vline(xintercept=c(x1min,x1max)) +
  geom_hline(yintercept=c(x2min,x2max)) +
  geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("GxE _extended") + xlab("z") + ylab("x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()


k1 <- length(x1grid)
k2 <- length(x2grid)
mu <- coef(obj1)$InterceptMB

eff1 <- B1x %*% theta1
eff2 <- B2x %*% theta2

pred1 <- data.frame(x1=x1grid, y=eff1)
pred2 <- data.frame(x2=x2grid, y=eff2)

ggplot(data=pred1, aes(x=x1, y=eff1)) + geom_line()
ggplot(data=pred2, aes(x=x2, y=eff2)) + geom_line()

yfit = kronecker(B1x,B2x) %*% theta3

x1 <- kronecker(x1grid, rep(1,k2))
x2 <- kronecker(rep(1,k1), x2grid)
pred <- data.frame(x1=x1,
                   x2=x2,
                   y = yfit)

ggplot(pred)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=y))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
  #  scale_fill_gradient(low = "blue", high = "red",name="Fitted") +
  geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("GxE") + xlab("z") + ylab("x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()

yfit = kronecker(rep(1,k1), eff2) + kronecker(B1x,B2x) %*% theta3
x1 <- kronecker(x1grid, rep(1,k2))
x2 <- kronecker(rep(1,k1), x2grid)
pred <- data.frame(x1=x1, x2=x2, y = yfit)

ggplot(pred, aes(x1,x2,z=y)) +
  geom_raster(mapping=aes(x=x1,y=x2,fill=y))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
  #  scale_fill_gradient(low = "blue", high = "red",name="Fitted") +
  geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("G + GxE") + xlab("z") + ylab("x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()

yfit = mu + kronecker(eff1,rep(1,k2)) + kronecker(rep(1,k1), eff2) + kronecker(B1x,B2x) %*% theta3

pred <- data.frame(x1=x1, x2=x2, y = yfit)

ggplot(pred)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=y))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
#  scale_fill_gradient(low = "blue", high = "red",name="Fitted") +
  geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("E + G + GxE") + xlab("z") + ylab("x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()

