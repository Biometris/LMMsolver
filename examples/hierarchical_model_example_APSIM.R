#
# Martin Boer, Biometris, Wageningen, the Netherlands.
#
# example of analysis of APSIM data using splines.
#
rm(list = ls())
library(LMMsolver)
library(dplyr)
library(ggplot2)

# subset of APSIM simulation data set Daniela Bustos-Korts
dat_traj = read.csv('APSIM_Emerald_1993_25genotypes.csv',header=TRUE, stringsAsFactors = FALSE)

head(dat_traj)

# use ton's/ha
dat_traj$biomass <- dat_traj$biomass/1000

# data for the simulations:
xmin =  30
xmax = 130
step =   5
#step =   1
dat <- filter(dat_traj, das %in% seq(xmin, xmax , by=step))
dim(dat)

#sel_geno <- paste0("g",formatC(c(8, 22),width=3,flag=0))
#dat <- filter(dat, geno %in% sel_geno)

set.seed(1234)
dat$z <- dat$das
dat$env <- as.factor(dat$das)
N <- nrow(dat)
sigma2e = 0.1
#sigma2e = 0.05
dat$ysim <- dat$biomass + rnorm(N,sd=sqrt(sigma2e))

labels <- unique(dat$geno)
dat$geno <- factor(dat$geno,levels=labels)
dat$env <- as.factor(dat$env)
dat$g_nr <- as.numeric(dat$geno)

Ngeno <- nlevels(dat$geno)
Nenv <- nlevels(dat$env)
Ngeno
Nenv

# how strong is the genetic signal?
obj0 <- LMMsolve(ysim ~ env, random=~geno, data=dat)
obj0$ED

# here we define the splines
# 1: E cubical splines, second order differences
# 2: G first degree splines, ridge penalty

## 1) definitions for environmenal covariate:

degr1 = 3
pord1 = 2
xmin1 <- xmin
xmax1 <- xmax
nseg1 <- 10
#nseg1 <- 25

knots1 = PsplinesKnots(xmin1, xmax1, degr1, nseg1)
B1 <- Bsplines(knots1, dat$z)
q1 <- ncol(B1)

U1sc <- calcUsc(q1, pord1)

## 1) definitions for genotype

# a very simple simple B-spline basis, helps that
# data can be in any order
knots2 = PsplinesKnots(1, Ngeno, 1, Ngeno-1)
B2 <- Bsplines(knots2, dat$g_nr)

# define matrix orthogonal to constant
J = matrix(data=1, nrow=Ngeno, ncol=Ngeno)
K = diag(Ngeno) - (1/Ngeno) * J
U2sc = eigen(K)$vectors[,-Ngeno, drop=FALSE]

# define the mixed model equations and solve:

lZ <- list()
lZ[[1]] = B1 %*% U1sc  #env
lZ[[2]] = B2 %*% U2sc  #geno
lZ[[3]] = dat$z*(B2 %*% U2sc) # x.geno
lZ[[4]] = RowKronecker(B1,B2) %*% kronecker(U1sc, U2sc)

Z <- do.call("cbind", lZ)
dat_ext = cbind(dat, Z)

lM <- ndxMatrix(dat, lZ, c("f(z)","g","g.z","f_g(z)"))
obj5 <- LMMsolve(ysim~z, randomMatrices=lM, data=dat_ext)
round(obj5$ED, 2)

z0 <- seq(xmin1,xmax1,by=1.0)
B1grid <- Bsplines(knots1, z0)
B1gridDz <- Bsplines(knots1, z0, deriv=TRUE)
mu <- coef(obj5)$'(Intercept)'
beta <- coef(obj5)$z
theta <- U1sc %*% coef(obj5)$'f(z)'
sum(theta)
ypredEnv   <- mu + beta*z0 + B1grid   %*% theta
ypredEnvDz <-      beta    + B1gridDz %*% theta
ypredE.df <- data.frame(z=z0, y=ypredEnv, dy=ypredEnvDz)
head(ypredE.df)

G_eff <- as.vector(U2sc %*% coef(obj5)$g)
G.z <- as.vector(U2sc %*% coef(obj5)$g.z)
ypredG <- as.vector(kronecker(rep(1,length(z0)), G_eff))
ypredGx <- as.vector(kronecker(z0, G.z))
ypredGxE <- ypredGx +
  kronecker(B1grid,diag(Ngeno)) %*% kronecker(U1sc,U2sc) %*% coef(obj5)$'f_g(z)'
ypredGtot <- ypredG + ypredGxE

M <- matrix(data=ypredGxE,nrow=length(z0),ncol=Ngeno,byrow=TRUE)
range(rowSums(M))

# calculate derivatives:
ypredGDz <- as.vector(kronecker(rep(1,length(z0)), G.z))
ypredGxEDz <- ypredGDz +
  kronecker(B1gridDz,diag(Ngeno)) %*% kronecker(U1sc,U2sc) %*% coef(obj5)$'f_g(z)'

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
p + geom_vline(xintercept=57)

p <- ggplot(dat, aes(x=z, y=ysim,col=geno)) + geom_point()
p + geom_line(ypredE.df, mapping=aes(x=z,y=y),col='black',size=1.5) +
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
p <- ggplot(dat_sel, aes(x=z, y=ysim,col=geno)) + geom_point() +
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


