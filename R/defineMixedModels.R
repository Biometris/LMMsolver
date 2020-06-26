makeZmatrix <- function(ran,dat)
{
  # random part of model, see implementation in SpATS....
  mf <- model.frame(ran, dat, drop.unused.levels = TRUE, na.action = NULL)
  mt <- terms(mf)    
  f.terms <- all.vars(mt)[attr(mt,"dataClasses") == "factor"]
  Z <- model.matrix(mt, data = mf, contrasts.arg = lapply(mf[,f.terms, drop = FALSE], contrasts, contrasts = FALSE))
  Z = Z[,-1, drop = FALSE]
  Z
}

makeGlist <- function(ran,dat)
{
  mf <- model.frame(ran, dat, drop.unused.levels = TRUE, na.action = NULL)
  mt <- terms(mf)    
  f.terms <- all.vars(mt)[attr(mt,"dataClasses") == "factor"]
  nlevelsRandom = sapply(mf,nlevels)
  
  # components of Ginv
  e <- cumsum(nlevelsRandom)
  s <- e - nlevelsRandom + 1
  lGinv = list()
  for (i in 1:length(s))
  {
    tmp <- rep(0,ncol(Z))
    tmp[(s[i]:e[i])] = 1.0
    lGinv[[i]] = cleanup(diag.spam(tmp))
  }
  names(lGinv) = names(mf)
  lGinv
}

# components of Rinv = phi1*A[[1]]+phi2*A[[2]]:
makeRlist <- function(df, column)
{
  lRinv = list()
  levels_f = levels(df[[column]])
  for (i in 1:length(levels_f)) {
    lRinv[[i]] = diag.spam(1.0*(df[[column]]==levels_f[i]))
  }
  names(lRinv) <- paste0(column,"_",levels_f,"!R")
  lRinv
}
