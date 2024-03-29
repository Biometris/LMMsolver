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
library(animation)
library(statgenSTA)
library(statgenGxE)

# load obj
solve_LMM = FALSE

# write to pdf or gif:
write_files = FALSE

#dat_all <- read.delim("SimulatedData124Envs_Biomass.txt")
#save(dat_all, file="SimulatedData124Envs_Biomass.rda")
load("SimulatedData124Envs_Biomass.rda")
dim(dat_all)

sel_env <- data.frame(env=c("Emerald_1993", "Emerald_2005",
                            "Merredin_2002","Merredin_2007",
                            "Narrabri_1993", "Narrabri_1998",
                            "Yanco_1993",	 "Yanco_1996",
                            "Narrabri_2003",	"Narrabri_2008",
                            "Narrabri_2011",	"Narrabri_2013"),
                      envtype = rep(LETTERS[1:3],each=4))
sel_env

dat_traj <- filter(dat_all, Env %in% sel_env$env) %>% droplevels()

Ngeno_sel <- 25
sel_geno <- paste0("g",formatC(c(1:Ngeno_sel),width=3,flag='0'))
dat_traj <- filter(dat_traj, geno %in% sel_geno) %>% droplevels()

# use ton's/ha
dat_traj$biomass <- dat_traj$biomass/1000

# data for the simulations:
xmin =  30
xmax = 160
step =   5
#step =   1
dat <- filter(dat_traj, das %in% seq(xmin, xmax , by=step))
dim(dat)

set.seed(1234)
N <- nrow(dat)
sigma2e = 0.1
#sigma2e = 0.05
dat$ysim <- dat$biomass + rnorm(N,sd=sqrt(sigma2e))

Glabels <- unique(dat$geno)
Elabels <- unique(dat$Env)
Glabels
Elabels

dat$g_nr <- as.numeric(dat$geno)
dat$e_nr <- as.numeric(dat$Env)

Ngeno <- nlevels(dat$geno)
Nenv  <- length(unique(dat$e_nr))
Ngeno
Nenv

# here we define the splines
# 1: time (days after sowing) cubical splines, second order differences
# 2: G first degree splines, ridge penalty
# 3: ENV (e_nr) first degree splines, ridge penalty

## 1) definitions for time (days after sowing)

degr1 = 3
pord1 = 2
xmin1 <- xmin
xmax1 <- xmax
nseg1 <- 12
#nseg1 <- 25

knots1 = PsplinesKnots(xmin1, xmax1, degr1, nseg1)
B1 <- as.spam(Bsplines(knots1, dat$das))
q1 <- ncol(B1)
q1

# a very simple simple B-spline basis, helps that
# data can be in any order
knots2 = PsplinesKnots(1, Ngeno, 1, Ngeno-1)
B2 <- as.spam(Bsplines(knots2, dat$g_nr))

knots3 = PsplinesKnots(1, Nenv, 1, Nenv-1)
B3 <- as.spam(Bsplines(knots3, dat$e_nr))

# linear space
# ~gitprojects/MBnotes/sparse_mixed_model_splines.tex/pdf
nknots1 <- length(knots1)
tau <- rollmean(knots1[-c(1,nknots1)], k=degr1)
# B1(z) %*% tau = z
all.equal(as.vector(B1 %*% tau), dat$das)

# intro example:
dat2 <- filter(dat, Env=='Narrabri_2003' & geno=='g001') %>%droplevels()
dat_traj2 <- filter(dat, Env=='Narrabri_2003' & geno=='g001') %>%droplevels()
D1 <- diff.spam(diag.spam(q1),    diff=pord1)


B1_2 <- as.spam(Bsplines(knots1, dat2$das))
BtB <- crossprod(B1_2)
DtD <- crossprod(D1)
BtY <- t(B1_2) %*% dat2$ysim
grid <- seq(4,-6, by=-0.5)
t0 <- seq(xmin1, xmax1, by=1.0)
B1grid <- as.spam(Bsplines(knots1, t0))

if (write_files) {
  movie.name = paste0("Psplines.gif")
  ani.options(loop = 1)
  saveGIF(
    for (i in 1:length(grid))
    {
      cat ("analysing penalty ", i, "\n")
      tmp <- grid[i]
      lambda = 10^tmp
      a <- solve(BtB+lambda*DtD,BtY)
      ypred <- B1grid %*% a
      df <- data.frame(das=t0,y=ypred)
      p <- ggplot(dat2, aes(x=das,y=ysim)) + geom_point(col='blue') +
        ggtitle(paste0("genotype g001 in environment Narrabri_2003, penalty = 10^",grid[i])) +
        xlab("days after sowing") + ylab("Biomass (ton/ha)") + ylim(0,15) +
        geom_line(dat_traj2,mapping=aes(x=das,y=biomass),col='green',size=1.0) +
        geom_line(df,mapping=aes(x=das,y=ypred),col='red',size=1.0) +
        theme(plot.title = element_text(hjust = 0.5))
      plot(p)
    },
    movie.name = movie.name)
}

lZ <- list()
lZ[[1]] = B1_2 %*% t(D1)  # das

Z <- as.matrix(do.call("cbind", lZ))
dat2_ext = cbind(dat2, Z)
lM <- ndxMatrix(dat2, lZ, c("f(t)"))

# define precision matrices:
DtD1 <- crossprod(D1)
precM1 <- D1 %*% DtD1 %*% t(D1)

lGinv <- list()
lGinv[['f(t)']]   <- precM1

obj0 <- LMMsolve(ysim~das, group=lM,lGinverse=lGinv, data=dat2_ext,tolerance=1.0e-4,
    display=TRUE,trace=TRUE)

t0 <- seq(xmin1, xmax1, by=1.0)
B1grid <- as.spam(Bsplines(knots1, t0))
B1gridDt <- as.spam(Bsplines(knots1, t0, deriv=TRUE))
mu <- coef(obj0)$'(Intercept)'
beta <- coef(obj0)$das
theta <- t(D1) %*% coef(obj0)$'f(t)'
sum(theta)
ypredMain   <- mu + beta*t0 + B1grid   %*% theta
ypredMainDt <-      beta    + B1gridDt %*% theta
ypred.df <- data.frame(das=t0,y=ypredMain)

# parameters for the B-spline basis:
a <- theta+beta*tau+mu
dfBbasis <- NULL
for (i in 1:q1)
{
  tmp <- data.frame(BsplineNr=i,t0=t0, B=B1grid[,i],Ba=B1grid[,i] * a[i],
                    dBdt = B1gridDt[,i], BaDt=B1gridDt[,i] * a[i])
  dfBbasis <- rbind(dfBbasis,tmp)
}
dfBbasis$BsplineNr <- as.factor(dfBbasis$BsplineNr)

p <- ggplot(dat2, aes(x=das,y=ysim)) + geom_point(col='blue') +
  ggtitle("Example: genotype g001 in environment Narrabri_2003") +
  xlab("days after sowing") + ylab("y") +
  geom_line(dfBbasis,mapping=aes(x=t0,y=Ba,col=BsplineNr)) +
  geom_line(ypred.df,mapping=aes(x=das,y=y)) +
  theme(plot.title = element_text(hjust = 0.5))
p


# sparse model:
D1 <- diff.spam(diag.spam(q1),    diff=pord1)
D2 <- diff.spam(diag.spam(Ngeno), diff=1)
D3 <- diff.spam(diag.spam(Nenv) , diff=1)

lZ <- list()
B1B2B3 <- RowKronecker(B1, RowKronecker(B2, B3))
lZ[[1]] = B1 %*% t(D1)  # das
lZ[[2]] = B2 %*% t(D2)  # geno
lZ[[3]] = B3 %*% t(D3)  # env
lZ[[4]] = RowKronecker(B1,B2) %*% kronecker(tau, t(D2))  # lin geno
lZ[[5]] = RowKronecker(B1,B3) %*% kronecker(tau, t(D3))  # lin env
lZ[[6]] = RowKronecker(B2,B3) %*% kronecker(t(D2),t(D3)) # GxE
lZ[[7]] = B1B2B3 %*% kronecker(tau, t(D2) %x% t(D3) )    # lin GxE
lZ[[8]] = RowKronecker(B1,B2) %*% kronecker(t(D1), t(D2)) # f_g(t)
lZ[[9]] = RowKronecker(B1,B3) %*% kronecker(t(D1), t(D3)) # f_g(t)
lZ[[10]] = B1B2B3 %*% (t(D1) %x% t(D2) %x% t(D3)) # f_g(t)

Z <- as.matrix(do.call("cbind", lZ))
dat_ext = cbind(dat, Z)
lM <- ndxMatrix(dat, lZ, c("f(t)","g","e","g.t","e.t","gxe","gxe.t","f_g(t)","f_e(t)",
                           "f_gxe(t)"))

# define precision matrices:
DtD1 <- crossprod(D1)
I_g <- diag.spam(1,Ngeno)
I_e <- diag.spam(1,Nenv)
precM1 <- D1 %*% DtD1 %*% t(D1)
precM2 <- D2 %*% I_g  %*% t(D2)
precM3 <- D3 %*% I_e  %*% t(D3)

lGinv <- list()
lGinv[['f(t)']]   <- precM1
lGinv[['g']]      <- precM2
lGinv[['e']]      <- precM3
lGinv[['g.t']]    <- precM2
lGinv[['e.t']]    <- precM3
lGinv[['gxe']]  <- kronecker(precM2, precM3)
lGinv[['gxe.t']]  <- kronecker(precM2, precM3)
lGinv[['f_g(t)']] <- kronecker(precM1,precM2)
lGinv[['f_e(t)']] <- kronecker(precM1,precM3)
lGinv[['f_gxe(t)']] <- precM1 %x% precM2 %x% precM3

if (solve_LMM)
{
  s <- proc.time()[3]
  obj <- LMMsolve(ysim~das, group=lM,lGinverse=lGinv, data=dat_ext,tolerance=1.0e-4,
                       display=TRUE,trace=TRUE)
  e <- proc.time()[3]
  cat("Computation time ", e-s, "seconds \n")
  save(obj, file="LMMsolve_APSIM_multi_env.rda")

} else {
  load(file="LMMsolve_APSIM_multi_env.rda")
}
round(obj$ED, 2)
obj$EDmax

EDtable0 <- data.frame(term=c('intercept','slope'), ED=c(1,1), EDmax=c(1,1), eff=c('F','F'))
EDtable <- data.frame(term=names(obj$ED), ED=round(obj$ED,2), EDmax=obj$EDmax, eff=rep('R',11))
EDtable <- rbind(EDtable0, EDtable)
EDtable <- EDtable[-3,]
rownames(EDtable) <- NULL
EDtableSum <- data.frame(term='TotalSum', ED=sum(EDtable$ED), EDmax=sum(EDtable$EDmax))
EDtable
EDtableSum

# make predictions on a dense grid:
t0 <- seq(xmin1, xmax1, by=1.0)
B1grid <- as.spam(Bsplines(knots1, t0))
B1gridDt <- as.spam(Bsplines(knots1, t0, deriv=TRUE))
mu <- coef(obj)$'(Intercept)'
beta <- coef(obj)$das
theta <- t(D1) %*% coef(obj)$'f(t)'
sum(theta)
ypredMain   <- mu + beta*t0 + B1grid   %*% theta
ypredMainDt <-      beta    + B1gridDt %*% theta
ypredMain.df <- data.frame(z=t0, y=ypredMain, dy=ypredMainDt)
head(ypredMain.df)

p <- ggplot(ypredMain.df, aes(x=z,y=y)) + geom_line(size=1.5,col='blue') +
  ggtitle("Mean Biomass as function of time") +
  xlab("days after sowing") + ylab("biomass") +
  theme(plot.title = element_text(hjust = 0.5))
p

# parameters for the B-spline basis:
a <- theta+beta*tau+mu
dfBbasis <- NULL
for (i in 1:q1)
{
  tmp <- data.frame(BsplineNr=i,t0=t0, B=B1grid[,i],Ba=B1grid[,i] * a[i],
                    dBdt = B1gridDt[,i], BaDt=B1gridDt[,i] * a[i])
  dfBbasis <- rbind(dfBbasis,tmp)
}
dfBbasis$BsplineNr <- as.factor(dfBbasis$BsplineNr)

dfBbasisSel <- filter(dfBbasis, BsplineNr %in% c(5,10))

p <- ggplot(dfBbasisSel, aes(x=t0,y=B,col=BsplineNr)) + geom_line() +
  ggtitle("example of two cubical B-splines") +
  xlab("days after sowing") + ylab("y") + ylim(0,1) +
  theme(plot.title = element_text(hjust = 0.5))
p

p <- ggplot(dfBbasis, aes(x=t0,y=B,col=BsplineNr)) + geom_line() +
  ggtitle("B-splines basis") +
  xlab("days after sowing") + ylab("y") + ylim(0,1) +
  theme(plot.title = element_text(hjust = 0.5)) + geom_hline(yintercept=1.0)
p


p <- ggplot(ypredMain.df, aes(x=z,y=y)) + geom_line(size=1.5,col='blue') +
  geom_line(dfBbasis, mapping=aes(x=t0,y=Ba,col=BsplineNr)) +
  ggtitle("Mean Biomass as function of time") +
  xlab("days after sowing") + ylab("biomass") +
  theme(plot.title = element_text(hjust = 0.5))
p


p <- ggplot(ypredMain.df, aes(x=z,y=dy)) + geom_line(size=1.5, col='blue') +
  ggtitle("Mean Growth rate") +
  xlab("days after sowing") + ylab("growth rate (ton/ha per day)") + ylim(0,0.1) +
  theme(plot.title = element_text(hjust = 0.5))
p


G_eff <- as.vector(t(D2) %*% coef(obj)$g)
G.t <- as.vector(t(D2) %*% coef(obj)$g.t)
ypredG <- as.vector(kronecker(rep(1,length(t0)), G_eff))
ypredGl<- as.vector(kronecker(t0, G.t))

# interaction term with sum to zero constraints:
theta_fgt <- kronecker(t(D1), t(D2)) %*% coef(obj)$'f_g(t)'
M <- matrix(data=theta_fgt,nrow=q1,ncol=Ngeno,byrow=TRUE)
range(rowSums(M))
range(colSums(M))

ypredGnl <- kronecker(B1grid,diag(Ngeno)) %*% theta_fgt
ypredGtot <- ypredG + ypredGl + ypredGnl

# calculate derivatives:
ypredGlDt <- as.vector(kronecker(rep(1,length(t0)), G.t))
ypredGDt <- ypredGlDt + kronecker(B1gridDt,diag(Ngeno)) %*% theta_fgt

geno <- rep(Glabels,times=length(t0))
x1_ext <- kronecker(t0, rep(1,Ngeno))
x2_ext <- kronecker(rep(1,length(t0)), c(1:Ngeno))
predG <- data.frame(geno,x1=x1_ext,x2=x2_ext, ypredGtot,
                   ypredTot = rep(ypredMain,each=Ngeno) + ypredGtot,
                   ypredGDt = ypredGDt,
                   ypredTotDt = rep(ypredMainDt,each=Ngeno) + ypredGDt)
head(predG)

mxBiomassAverage <- max(ypredMain.df$y)
mxBiomass <- predG %>% group_by(geno) %>% summarize(mx=max(ypredTot))
mxBiomass$type <- as.factor(ifelse(mxBiomass$mx>mxBiomassAverage,"High","Low"))
predG <- left_join(predG, mxBiomass,by='geno')
head(predG)

dat <- left_join(dat,sel_env,by=c('Env'='env'))

p <- ggplot(dat, aes(x=das,y=ysim,col=envtype)) +
  facet_wrap(~geno) + geom_point() +
  ggtitle("Raw data for 12 environments") +
  xlab("days after sowing") + ylab("biomass") +
  theme(plot.title = element_text(hjust = 0.5))
p

p <- ggplot(predG, aes(x=x1, y=ypredTot,col=type)) +
  facet_wrap(~geno) + geom_line() +
  geom_line(ypredMain.df, mapping=aes(x=z,y=y),col='black',linetype='dashed') +
  ggtitle("Average biomass across 12 environments") + xlab("days after sowing") + ylab("biomass") +
  theme(plot.title = element_text(hjust = 0.5))
p

p <- ggplot(predG, aes(x=x1, y=ypredTotDt,col=type)) +
  facet_wrap(~geno) + geom_line() +
  ggtitle("Average growth rates across 12 environments") +
  geom_line(ypredMain.df, mapping=aes(x=z,y=dy),col='black',linetype='dashed') +
    xlab("days after sowing") + ylab("growth rate") +
  theme(plot.title = element_text(hjust = 0.5))
p

E_eff <- as.vector(t(D3) %*% coef(obj)$e)
E.t <- as.vector(t(D3) %*% coef(obj)$e.t)
ypredE <- as.vector(kronecker(rep(1,length(t0)), E_eff))
ypredEl<- as.vector(kronecker(t0, E.t))

theta_fEt <- kronecker(t(D1), t(D3)) %*% coef(obj)$'f_e(t)'
ypredEnl <- kronecker(B1grid,diag(Nenv)) %*% theta_fEt
ypredEtot <- ypredE + ypredEl + ypredEnl

# calculate derivatives:
ypredElDt <- as.vector(kronecker(rep(1,length(t0)), E.t))
ypredEDt <- ypredElDt + kronecker(B1gridDt,diag(Nenv)) %*% theta_fEt

env <- rep(Elabels,times=length(t0))
x1_ext <- kronecker(t0, rep(1,Nenv))
x2_ext <- kronecker(rep(1,length(t0)), c(1:Nenv))
predE <- data.frame(env,x1=x1_ext,x2=x2_ext,
                   ypredTot = rep(ypredMain,each=Nenv) + ypredEtot,
                   ypredTotDt = rep(ypredMainDt,each=Nenv) + ypredEDt)
head(predE)

predE <- left_join(predE,sel_env,by='env')

p <- ggplot(predE, aes(x=x1, y=ypredTot,col=envtype)) +
  facet_wrap(~env) + geom_line() +
  geom_line(ypredMain.df, mapping=aes(x=z,y=y),col='black',linetype='dashed') +
  ggtitle("Average biomass of 25 genotypes") + xlab("days after sowing") + ylab("biomass") +
  theme(plot.title = element_text(hjust = 0.5))
p

p <- ggplot(predE, aes(x=x1, y=ypredTotDt,col=envtype)) +
  facet_wrap(~env) + geom_line() +
  geom_line(ypredMain.df, mapping=aes(x=z,y=dy),col='black',linetype='dashed') +
  ggtitle("Average growth rates of 25 genotypes") + xlab("days after sowing") + ylab("growth rate") +
  theme(plot.title = element_text(hjust = 0.5))
p

# make full time x G x E predictions

# vector of ones..
oneT <- rep(1,q1)
oneE <- rep(1,Nenv)
oneG <- rep(1,Ngeno)

Glabels <- as.character(Glabels)
Elabels <- as.character(Elabels)

x1 = rep(t0, each=Ngeno*Nenv)
x2 = rep(rep(Glabels,each=Nenv),times=length(t0))
x3 = rep(Elabels,times=length(t0)*Ngeno)

B1B2B3   <- as.spam(B1grid %x% diag(1,Ngeno) %x% diag(1,Nenv))
B1B2B3Dt <- as.spam(B1gridDt %x% diag(1,Ngeno) %x% diag(1,Nenv))


lM <- ndxMatrix(dat, lZ, c("f(t)","g","e","g.t","e.t","gxe","gxe.t","f_g(t)","f_e(t)",
                           "f_gxe(t)"))

D3<-as.matrix(D3)
pred0 <- mu + beta*x1
pred1 <- B1B2B3 %*% (t(D1) %x% oneG %x% oneE) %*% coef(obj)$'f(t)'
pred2 <- B1B2B3 %*% (oneT  %x% t(D2) %x% oneE) %*% coef(obj)$g
pred3 <- B1B2B3 %*% (oneT  %x% oneG %x% t(D3)) %*% coef(obj)$e

pred4 <- B1B2B3 %*% (tau  %x% t(D2) %x% oneE) %*% coef(obj)$'g.t'
pred5 <- B1B2B3 %*% (tau  %x% oneG %x% t(D3)) %*% coef(obj)$'e.t'

pred6 <- B1B2B3 %*% (oneT  %x% t(D2) %x% t(D3)) %*% coef(obj)$'gxe'
pred7 <- B1B2B3 %*% (tau  %x% t(D2) %x% t(D3)) %*% coef(obj)$'gxe.t'

pred8 <- B1B2B3 %*% (t(D1)  %x% t(D2) %x% oneE) %*% coef(obj)$'f_g(t)'
pred9 <- B1B2B3 %*% (t(D1)  %x% oneG %x% t(D3)) %*% coef(obj)$'f_e(t)'
pred10 <- B1B2B3 %*% (t(D1)  %x% t(D2) %x% t(D3)) %*% coef(obj)$'f_gxe(t)'

ypred <- pred0+pred1+pred2+pred3+pred4+pred5+pred6+pred7+pred8+pred9+pred10

pred0 <- beta
pred1 <- B1B2B3Dt %*% (t(D1) %x% oneG %x% oneE) %*% coef(obj)$'f(t)'
pred2 <- B1B2B3Dt %*% (oneT  %x% t(D2) %x% oneE) %*% coef(obj)$g
pred3 <- B1B2B3Dt %*% (oneT  %x% oneG %x% t(D3)) %*% coef(obj)$e

pred4 <- B1B2B3Dt %*% (tau  %x% t(D2) %x% oneE) %*% coef(obj)$'g.t'
pred5 <- B1B2B3Dt %*% (tau  %x% oneG %x% t(D3)) %*% coef(obj)$'e.t'

pred6 <- B1B2B3Dt %*% (oneT  %x% t(D2) %x% t(D3)) %*% coef(obj)$'gxe'
pred7 <- B1B2B3Dt %*% (tau  %x% t(D2) %x% t(D3)) %*% coef(obj)$'gxe.t'

pred8 <- B1B2B3Dt %*% (t(D1)  %x% t(D2) %x% oneE) %*% coef(obj)$'f_g(t)'
pred9 <- B1B2B3Dt %*% (t(D1)  %x% oneG %x% t(D3)) %*% coef(obj)$'f_e(t)'
pred10 <- B1B2B3Dt %*% (t(D1)  %x% t(D2) %x% t(D3)) %*% coef(obj)$'f_gxe(t)'

ypredDt <- pred0+pred1+pred2+pred3+pred4+pred5+pred6+pred7+pred8+pred9+pred10

pred <- data.frame(das=x1,geno=x2,env=x3, ypred=ypred,ypredDt=ypredDt,
                   ypredGxEDt=pred6+pred7+pred10)

#sel_geno <- paste0("g",formatC(c(9,12,17),width=3,flag=0))
sel_geno <- paste0("g",formatC(c(1:25),width=3,flag=0))
#sel_geno <- paste0("g",formatC(c(7,8),width=3,flag=0))

pred_sel <- filter(pred, geno %in% sel_geno)


p <- ggplot(pred_sel, aes(x=das, y=ypred,col=geno)) +
  facet_wrap(~env) + geom_line() +
  ggtitle("Biomass") + xlab("days after sowing") + ylab("Biomass") +
  theme(plot.title = element_text(hjust = 0.5))
p


p <- ggplot(pred_sel, aes(x=das, y=ypredDt,col=geno)) +
  facet_wrap(~env) + geom_line() +
  ggtitle("Growth rate") + xlab("days after sowing") + ylab("growth rate") +
  theme(plot.title = element_text(hjust = 0.5))
p

#p <- ggplot(pred_sel, aes(x=das, y=ypredDt,col=env)) +
#  facet_wrap(~geno) + geom_line() +
#  ggtitle("Average biomass of 25 genotypes") + xlab("days after sowing") + ylab("biomass") +
#  theme(plot.title = element_text(hjust = 0.5))
#p

#dat_traj_sel <- filter(dat_traj,geno %in% sel_geno)

#p <- ggplot(dat_traj_sel, aes(x=das, y=biomass,col=geno)) +
#  facet_wrap(~Env) + geom_line() +
#  ggtitle("Average biomass of 25 genotypes") + xlab("days after sowing") + ylab("biomass") +
#  theme(plot.title = element_text(hjust = 0.5))
#p


sel_day <- 90
pred_sel_day <- filter(pred,das==sel_day)

pred_sel_day <- left_join(pred_sel_day,sel_env,by='env')
pred_sel_day$genoGrp <- as.factor(ifelse(pred_sel_day$geno=='g008','no','yes'))

TDobj <- statgenSTA::createTD(data=pred_sel_day,genotype='geno', trial='env')
AMMIobj <- gxeAmmi(TD = TDobj, trait = "ypredDt")

pred_sel <- filter(pred_sel, geno=='g008')
pred_sel <- left_join(pred_sel,sel_env,by='env')

p <- ggplot(pred_sel, aes(x=das, y=ypredGxEDt,col=envtype)) +
  facet_wrap(~env) + geom_line(size=1.0) +
  ggtitle("GxE for Growth rate genotype g008") + xlab("days after sowing") + ylab("growth rate") +
  geom_vline(xintercept=90,linetype='dashed') + geom_vline(xintercept=120,linetype='dashed') +
  theme(plot.title = element_text(hjust = 0.5))
p


p <- plot(AMMIobj, scale=0.5, plotType="AMMI2", sizeGeno=3, colorEnvBy='envtype',
   colorGenoBy='genoGrp', colEnv=c('red','green','blue'), colGeno=c('blue','black'),
    output=FALSE)
p +  ggtitle(paste0("AMMI 2 Growth rate at day ", sel_day)) + theme(legend.position = "none")

#colorGenoBy='GenoGrp',colorEnvBy = 'EnvGrp',
#colEnv = c('blue','green', 'red'),
#colGeno=c('blue','green','black'))


p <- ggplot(pred_sel, aes(x=das, y=ypredDt,col=geno)) +
  facet_wrap(~env) + geom_line() +
  ggtitle("Growth rate") + xlab("days after sowing") + ylab("growth rate") +
  geom_vline(xintercept=90,linetype='dashed') + geom_vline(xintercept=120,linetype='dashed') +
  theme(plot.title = element_text(hjust = 0.5))
p

if (write_files) {
  movie.name = paste0("AMMI.gif")
  ani.options(loop = 1)
  saveGIF(
    for (i in 90:120)
    {
      cat ("analysing day ", i, "\n")
      predday <- filter(pred, das==i)
      predday <- left_join(predday,sel_env,by='env')
      predday$genoGrp <- as.factor(ifelse(predday$geno=='g008','no','yes'))

      TDobj <- statgenSTA::createTD(data=predday,genotype='geno', trial='env')
      AMMIobj <- gxeAmmi(TD = TDobj, trait = "ypredDt")

      p <- plot(AMMIobj, scale=0.5, plotType="AMMI2", sizeGeno=3, colorEnvBy='envtype',
                colorGenoBy='genoGrp', colEnv=c('red','green','blue'), colGeno=c('blue','black'),
                output=FALSE)
      p <- p +  ggtitle(paste0("AMMI 2 Growth rate at day ", i)) + theme(legend.position = "none")
      plot(p)
    },
    movie.name = movie.name)
}






