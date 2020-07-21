#' ---
#' title: GxE analysing using ANOVA smoothing
#' author: Martin Boer
#' ---

rm(list=ls())

library(ggplot2)
library(sampling)
library(LMMsolver)

set.seed(1234)

# residual error
sigma2 <- 25.0

# simulated surface a rotated ellipse
simfun <- function(x1,x2, theta) {
  a = 0.3
  b = 0.2
  x1c = 0.4
  x2c = 0.6
  x1r <- cos(theta)*(x1-x1c) - sin(theta)*(x2-x2c)
  x2r <- sin(theta)*(x1-x1c) + cos(theta)*(x2-x2c)
  25 - (x1r^2/a^2 + x2r^2/b^2)
}

# clockwise rotation of the simulated ellipse surface:
degrees_rotation = -30
rotation = (degrees_rotation/180)*pi

x1grid <- (c(1:100)-0.5)/100
x2grid <- (c(1:100)-0.5)/100
grid <- expand.grid(x1grid, x2grid)
names(grid) <- c("x1","x2")

# 20 environments:
# 25 genotypes:
nEnv = 20
nGeno = 25

x1f <- findInterval(x1grid,seq(from=0,to=1,length=nEnv+1))
ndx1 <- strata(data.frame(x1f),stratanames="x1f",size=rep(1,times=nEnv),method='srswor')$ID_unit

x2f <- findInterval(x1grid,seq(from=0,to=1,length=nGeno+1))
ndx2 <- strata(data.frame(x2f),stratanames="x2f",size=rep(1,times=nGeno),method='srswor')$ID_unit

grid$ysim <- simfun(grid$x1, grid$x2, theta=rotation)

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
sim.df$ysim <- simfun(sim.df$x1, sim.df$x2, theta= rotation)
e <- rnorm(Nobs, mean = 0, sd = sqrt(sigma2))
sim.df$y <- sim.df$ysim + e

# add env and genotype number:
sim.df$geno=rep(1:nGeno,each=nEnv)
sim.df$env=rep(1:nEnv,times=nGeno)

ggplot(sim.df)+
  geom_raster(mapping=aes(x=env,y=geno,fill=y))+
  scale_fill_gradientn(name="ysim",colours= topo.colors(100)) +
  #geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("simulated") + xlab("env") + ylab("genotype") +
  #scale_x_continuous(breaks=c(1,6,nEnv)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()


degree = 3
pord = 1
x1min = 0
x1max = 1
x2min = 0
x2max = 1
nseg1 = 10
nseg2 = 8

knots1 = PsplinesKnots(x1min, x1max, degree, nseg1)
knots2 = PsplinesKnots(x2min, x2max, degree, nseg2)


x1 = sim.df$x1
x2 = sim.df$x2

B1  <- Bsplines(knots1, x1)
B2  <- Bsplines(knots2, x2)
B12 <- RowKronecker(B1,B2)

q1 <- ncol(B1)
q2 <- ncol(B2)

U1sc <- calcUsc(q1, pord)
U2sc <- calcUsc(q2, pord)
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

# model without interaction:

lZ <- list()
lZ[[1]] = B1 %*% U1sc
lZ[[2]] = B2 %*% U2sc
Z <- do.call("cbind", lZ)

df_ext = cbind(df, Z)

lM <- ndxMatrix(df, lZ, c("E","G"))
obj0 = LMMsolve(fixed=y~1, random=NULL,randomMatrices=lM, data=df_ext)
obj0$logL
obj0$ED

# model with GxE interaction:
lZ <- list()
lZ[[1]] = B1  %*% U1sc
lZ[[2]] = B2  %*% U2sc
lZ[[3]] = B12 %*% U12sc
Z <- do.call("cbind",lZ)
df_ext = cbind(df,Z)

lM <- ndxMatrix(df, lZ, c("E","G","GxE"))
obj1 = LMMsolve(fixed=y~1, random=NULL,randomMatrices=lM, data=df_ext)
obj1$logL
obj1$ED

# just check the sum to zero constraints:

theta1 <- U1sc %*% coef(obj1)$E
theta2 <- U2sc %*% coef(obj1)$G
sum(theta1)
sum(theta2)
theta3 <- U12sc %*% coef(obj1)$GxE

# marginal equal to zero:
M = matrix(data=theta3,nrow=q1,ncol=q2,byrow=TRUE)
rowSums(M)
colSums(M)

# make predictions on a dense grid:
B1x = Bsplines(knots1, x1grid)
B2x = Bsplines(knots2, x2grid)

k1 <- length(x1grid)
k2 <- length(x2grid)
mu <- coef(obj1)$'(Intercept)'

eff1 <- B1x %*% theta1
eff2 <- B2x %*% theta2

pred1 <- data.frame(x1=x1grid, y=eff1)
pred2 <- data.frame(x2=x2grid, y=eff2)

ggplot(data=pred1, aes(x=x1, y=eff1)) + geom_line()
ggplot(data=pred2, aes(x=x2, y=eff2)) + geom_line()

# calculate the predictions on the grid:
fit_E <- kronecker(eff1,rep(1,k2))
fit_G <- kronecker(rep(1,k1), eff2)
fit_GxE = kronecker(B1x,B2x) %*% theta3
fit_GGE = fit_G + fit_GxE
yfit = mu + fit_E + fit_G + fit_GxE

x1_ext <- kronecker(x1grid, rep(1,k2))
x2_ext <- kronecker(rep(1,k1), x2grid)
pred <- data.frame(x1=x1_ext,x2=x2_ext, fit_G, fit_GxE, fit_GGE, yfit)

ggplot(pred)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=fit_GxE))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
  geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("GxE") + xlab("z") + ylab("x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()

ggplot(pred) +
  geom_raster(mapping=aes(x=x1,y=x2,fill=fit_GGE))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
  geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("G + GxE") + xlab("z") + ylab("x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()

ggplot(pred)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=yfit))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
  geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("E + G + GxE") + xlab("z") + ylab("x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()

ggplot(grid)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=ysim))+
  scale_fill_gradientn(name="ysim",colours= topo.colors(100)) +
  #  scale_fill_gradient(low = "blue", high = "red",name="z") +
  #geom_vline(xintercept=s1bnd) +
  #geom_hline(yintercept=s1bnd) +
  geom_point(sim.df,mapping=aes(x=x1,y=x2),size=2)+
  ggtitle("true simulated surface") + xlab("z") + ylab("x") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_equal()
