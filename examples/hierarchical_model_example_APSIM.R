#
# Martin Boer, Biometris, Wageningen, the Netherlands.
#
# example of analysis of APSIM data using splines.
#
rm(list = ls())
library(LMMsolver)
library(dplyr)
library(zoo)
library(ggplot2)

# subset of APSIM simulation data set Daniela Bustos-Korts
dat_traj = read.csv('APSIM_Emerald.csv',header=TRUE, stringsAsFactors = FALSE)

head(dat_traj)

Ngeno_sel <- 25
Year_sel <- 1993

sel_geno <- paste0("g",formatC(c(1:Ngeno_sel),width=3,flag='0'))
dat_traj <- filter(dat_traj, geno %in% sel_geno & year == Year_sel)

# use ton's/ha
dat_traj$biomass <- dat_traj$biomass/1000

# data for the simulations:
xmin =  30
xmax = 130
step =   5
#step =   1
dat <- filter(dat_traj, das %in% seq(xmin, xmax , by=step))
dim(dat)

set.seed(1234)
N <- nrow(dat)
sigma2e = 0.1
#sigma2e = 0.05
dat$ysim <- dat$biomass + rnorm(N,sd=sqrt(sigma2e))

labels <- unique(dat$geno)
dat$geno <- factor(dat$geno,levels=labels)
dat$g_nr <- as.numeric(dat$geno)

Ngeno <- nlevels(dat$geno)
Ngeno

# here we define the splines
# 1: time (days after sowing) cubical splines, second order differences
# 2: G first degree splines, ridge penalty

## 1) definitions for time (days after sowing)

degr1 = 3
pord1 = 2
xmin1 <- xmin
xmax1 <- xmax
nseg1 <- 10
#nseg1 <- 25

knots1 = PsplinesKnots(xmin1, xmax1, degr1, nseg1)
B1 <- as.spam(Bsplines(knots1, dat$das))
q1 <- ncol(B1)

# a very simple simple B-spline basis, helps that
# data can be in any order
knots2 = PsplinesKnots(1, Ngeno, 1, Ngeno-1)
B2 <- as.spam(Bsplines(knots2, dat$g_nr))

# linear space
# ~gitprojects/MBnotes/sparse_mixed_model_splines.tex/pdf
nknots1 <- length(knots1)
tau <- rollmean(knots1[-c(1,nknots1)], k=degr1)
# B1(z) %*% tau = z
all.equal(as.vector(B1 %*% tau), dat$das)

# sparse model:
D1 <- diff.spam(diag.spam(q1),    diff=pord1)
D2 <- diff.spam(diag.spam(Ngeno), diff=1)

lZ <- list()
lZ[[1]] = B1 %*% t(D1)  #env
lZ[[2]] = B2 %*% t(D2)  #geno
lZ[[3]] = RowKronecker(B1,B2) %*% kronecker(tau, t(D2))
lZ[[4]] = RowKronecker(B1,B2) %*% kronecker(t(D1), t(D2))

Z <- as.matrix(do.call("cbind", lZ))
dat_ext = cbind(dat, Z)

lM <- ndxMatrix(dat, lZ, c("f(t)","g","g.t","f_g(t)"))
I_g <- diag.spam(1,Ngeno)
DtD1 <- crossprod(D1)
lGinv <- list()
lGinv[['f(t)']] <- D1 %*% DtD1 %*% t(D1)
lGinv[['g']]   <- D2 %*% I_g %*% t(D2)
lGinv[['g.t']] <- D2 %*% I_g %*% t(D2)
lGinv[['f_g(t)']] <- kronecker(D1 %*% DtD1 %*% t(D1), D2 %*% I_g %*% t(D2))
names(lGinv)
obj <- LMMsolve(ysim~das, randomMatrices=lM,lGinverse=lGinv, data=dat_ext,
                       display=TRUE,monitor=TRUE)
round(obj$ED, 2)

z0 <- seq(xmin1,xmax1,by=1.0)
B1grid <- as.spam(Bsplines(knots1, z0))
B1gridDz <- as.spam(Bsplines(knots1, z0, deriv=TRUE))
mu <- coef(obj)$'(Intercept)'
beta <- coef(obj)$das
theta <- t(D1) %*% coef(obj)$'f(t)'
sum(theta)
ypredEnv   <- mu + beta*z0 + B1grid   %*% theta
ypredEnvDz <-      beta    + B1gridDz %*% theta
ypredE.df <- data.frame(z=z0, y=ypredEnv, dy=ypredEnvDz)
head(ypredE.df)

G_eff <- as.vector(t(D2) %*% coef(obj)$g)
G.t <- as.vector(t(D2) %*% coef(obj)$g.t)
ypredG <- as.vector(kronecker(rep(1,length(z0)), G_eff))
ypredGx <- as.vector(kronecker(z0, G.t))
ypredGxE <- ypredGx +
  kronecker(B1grid,diag(Ngeno)) %*% kronecker(t(D1), t(D2)) %*% coef(obj)$'f_g(t)'
ypredGtot <- ypredG + ypredGxE

M <- matrix(data=ypredGxE,nrow=length(z0),ncol=Ngeno,byrow=TRUE)
range(rowSums(M))

# calculate derivatives:
ypredGDz <- as.vector(kronecker(rep(1,length(z0)), G.t))
ypredGxEDz <- ypredGDz +
  kronecker(B1gridDz,diag(Ngeno)) %*% kronecker(t(D1),t(D2)) %*% coef(obj)$'f_g(t)'

geno <- rep(labels,times=length(z0))
x1_ext <- kronecker(z0, rep(1,Ngeno))
x2_ext <- kronecker(rep(1,length(z0)), c(1:Ngeno))
pred <- data.frame(geno,x1=x1_ext,x2=x2_ext, ypredGtot,
                   ypredTot = rep(ypredEnv,each=Ngeno) + ypredGtot,
                   ypredGxEDz = ypredGxEDz,
                   ypredTotDz = rep(ypredEnvDz,each=Ngeno) + ypredGxEDz)
head(pred)

# some plots...

p <- ggplot(dat_traj, aes(x=das, y=biomass,col=geno)) + geom_line() +
  ggtitle("APSIM deterministic simulations") + xlab("days of sowing") + ylab("biomass") +
  theme(plot.title = element_text(hjust = 0.5))
p

p <- ggplot(dat, aes(x=das, y=ysim,col=geno)) + geom_point()
p + geom_line(ypredE.df, mapping=aes(x=z,y=y),col='black',size=1.5) +
  xlab("days after sowing") +
  ggtitle("population mean curve and raw data") +  theme(plot.title = element_text(hjust = 0.5))

p <- ggplot(ypredE.df, aes(x=z, y=dy)) + geom_line() +
  ggtitle("mean growth rate") +  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("days after sowing") + ylab("growth rate [ton/ha per day]")
p

ggplot(pred)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=ypredGtot))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
  ggtitle("Genetic effects, deviations from mean") + xlab("days after sowing") +
  ylab("genotype") + theme(plot.title = element_text(hjust = 0.5))


ggplot(pred)+
  geom_raster(mapping=aes(x=x1,y=x2,fill=ypredGxEDz))+
  scale_fill_gradientn(name="Fitted",colours=topo.colors(100)) +
  ggtitle("Growth rates") + xlab("days after sowing") +
  ylab("genotype") + theme(plot.title = element_text(hjust = 0.5))

# selection of genotypes:
# genotype 6 has low biomass
# genotype 8 is late mature, still good biomass at the end
# genotype 22 is the best genotype in this environment
sel_geno <- paste0("g",formatC(c(6, 8, 22),width=3,flag=0))
pred_sel <- filter(pred, geno %in% sel_geno)

dat_sel <- filter(dat, geno %in% sel_geno)
dat_traj_sel <- filter(dat,geno %in% sel_geno)
p <- ggplot(dat_sel, aes(x=das, y=ysim,col=geno)) + geom_point() +
  #geom_line(ypredE.df, mapping=aes(x=z,y=y),col='black') +
  geom_line(pred_sel,mapping=aes(x=x1,y=ypredTot,col=geno)) +
  geom_line(dat_traj_sel, mapping=aes(x=das,y=biomass,col=geno), linetype='dashed') +
  ggtitle("Trajectories for selection of genotypes (dashed true curves)") +
  theme(plot.title = element_text(hjust = 0.5))
p

p <-ggplot(ypredE.df, aes(x=z, y=dy)) + geom_line() +
    geom_line(pred_sel, mapping=aes(x=x1,y=ypredTotDz,col=geno)) +
    ggtitle("growth rate") +  theme(plot.title = element_text(hjust = 0.5)) +
    xlab("days after sowing") + ylab("growth rate [ton/ha per day]")
p


