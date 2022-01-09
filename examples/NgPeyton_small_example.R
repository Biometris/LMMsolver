#' ---
#' title: Sparse Cholesky Decomposition
#' author: Martin Boer, Biometris
#' ---
#
# see paper Furrer and Sain, sept 2010,
# with correction in Figure 1.
#
rm(list = ls())
library(spam)
library(LMMsolver)

split = function(x)
{
  s = x[c(1:(length(x)-1))]
  e = x[-1] - 1
  df = data.frame(s,e)
  df
}

get.row = function(A,r)
{
  s = A@rowpointers[r]
  e = A@rowpointers[r+1] - 1
  range = c(s:e)
  df = data.frame(row = r, col=A@colindices[range], values = A@entries[range])
  df
}

# just for illustration, how we can easily get elements from A, a spam matrix.
convert.matrix = function(A)
{
  row_range = split(A@rowpointers)
  dim = A@dimension
  nrow = dim[1]
  ncol = dim[2]
  B = matrix(0,nrow=nrow,ncol=ncol)
  for (r in 1:nrow)
  {
    range = c(row_range$s[r]:row_range$e[r])
    col_nr = A@colindices[range]
    values = A@entries[range]
    B[r,col_nr] = values
  }
  B
}

# just for illustration, how we can get elements from U, a spam.chol.NgPeyton matrix
# print as in Ng and Peyton paper
convert.cholesky = function(U)
{
  row_range = split(U@rowpointers)
  col_range = split(U@colpointers)
  sno_range = split(U@supernodes)

  dim = U@dimension
  B = matrix(0,nrow=dim[1],ncol=dim[2])

  nsupernodes = nrow(sno_range)
  for (node in 1:nsupernodes)
  {
    cat("Node ", node, "\n" )
    rows.node = c(sno_range$s[node]:sno_range$e[node])
    cols.node = c(col_range$s[node]:col_range$e[node])
    for (r in rows.node)
    {
      cat(" col", r, "\n")
      cat(" rows", cols.node, "\n")

      row.range = c(row_range$s[r]:row_range$e[r])
      col_nr = U@colindices[cols.node]
      cat(" rowindices", col_nr, "\n")

      values = U@entries[row.range]
      B[r,col_nr] = values
      cols.node = cols.node[-1]
    }
  }
  B
}

# logdeterminant, for given Cholesky U:
log.det = function(U) {
  ndx = U@rowpointers[-length(U@rowpointers)]
  logdet = 0.0
  for (i in ndx) {
    logdet = logdet + 2.0*log(U@entries[i])
  }
  logdet
}

# example, see Figure 1 in Furrer and Sain paper (with correction):
A = matrix(c(1,   0.5,   0, 0.5, 0.5,
             0.5,   1,   0, 0.5, 0,
             0,     0,   1,   0, 0.5,
             0.5, 0.5,   0,   1, 0,
             0.5,   0, 0.5,   0, 1),ncol=5,byrow=TRUE)
A
A = as.spam(A)
display(A)
slotNames(A)
length(A@entries)
A@colindices
A@rowpointers

nelem.row = diff(A@rowpointers)
nelem.row

# sparse Cholesky Decomposition:
U = chol(A)
display(t(U))
slotNames(U)

# as shown in paper:
U@colindices
U@colpointers
U@rowpointers

# print supernodes:
U@supernodes

# get row 2 and 3 as examples:
get.row(A,2)
get.row(A,3)

# test conversion function:
B = convert.matrix(A)
all.equal(B, as.matrix(A))

# test conversion function:
C = convert.cholesky(U)
all.equal(C,as.matrix(U))

# nice small example with four supernodes:
n1 = 2
n2 = 5
D1 = diff(diag(n1),diff=1)
D2 = diff(diag(n2),diff=1)

A1 = diag(n1) + t(D1)%*%D1
A2 = diag(n2) + t(D2)%*%D2
B = kronecker(A1,A2)
B
B = as.spam(B)
U = chol(B)
U@supernodes
diff(U@supernodes)
display(t(U))

U@pivot
U@pivot[U@invpivot]
U@invpivot[U@pivot]

# take U and reorder back to B:
U2 = as.spam(U)
U2[U@pivot,U@pivot] = U2
B_rec = t(U2) %*% U2
all.equal(B,B_rec)

# reorder B, and compare with t(U) %*% U
B_ord = B[U@pivot,U@pivot]
U3 = as.spam(U)
B_rec2 = t(U3) %*% U3
all.equal(B_ord,B_rec2)

as.double(determinant(B)$modulus)
log.det(U)

# row i in U
i = 4
display(U)
p = U@pivot
pinv = U@invpivot
df = get.row(B,p[i])
df$ord = pinv[df$col]
df = df[order(df$ord), ]
df = df[df$ord>=i,]
df

U@supernodes
U@colpointers
U@rowpointers
U@colindices
U@entries

row_range = split(U@rowpointers)
col_range = split(U@colpointers)
sno_range = split(U@supernodes)
row_range
col_range
sno_range

display(t(U))
convert.cholesky(U)

lP <- list()
lP[[1]] <- B
ADchol <- LMMsolver:::ADchol(lP)
slotNames(ADchol)
ADchol@colpointers
ADchol@rowindices

ADcholnew <- LMMsolver:::ADcholnew(lP)

slotNames(ADcholnew)

ADcholnew@supernodes
ADcholnew@rowpointers
ADcholnew@colpointers
ADcholnew@rowindices

LMMsolver:::PrintADchol(ADcholnew, lambda=1.0)

display(t(U))
