// Backwards Automatic Differentiation of the Cholesky Algorithm
// for calculating partial derivatives of the log-determinant of
// positive definite symmetric sparse matrices.
//
// Implementation by Martin Boer, 2026.
//
// The algorithm combines the sparse Cholesky factorization of
// Ng and Peyton (1993) with the reverse differentiation approach
// of Smith (1995).
//
// References:
//
// Ng, E. G. and B. W. Peyton (1993).
// "Block sparse Cholesky algorithms on advanced uniprocessor computers."
// SIAM Journal on Scientific Computing 14, 1034-1056.
//
// Furrer, R. and S. R. Sain (2010).
// "spam: A sparse matrix R package with emphasis on MCMC
// methods for Gaussian Markov random fields."
// Journal of Statistical Software 36, 1-25.
//
// Smith, S. P. (1995).
// "Differentiation of the Cholesky algorithm."
// Journal of Computational and Graphical Statistics 4, 134-147.
//
// Smith, S. P. (2000).
// "A Tutorial on Simplicity and Computational Differentiation
// for Statisticians."

#include <Rcpp.h>
#include <set>
#include <vector>
#include "AuxFun.h"
#include "SparseMatrix.h"
#include "cholesky.h"
#include "ADcholesky.h"

using namespace Rcpp;
using namespace std;

// j is current column in Supernode J
void ADcmod1(NumericVector& F,
             const NumericVector& L, int j, int J,
             const IntegerVector& supernodes,
             const IntegerVector& colpointers)
{
  int s = colpointers[j];
  int e = colpointers[j+1];
  // for all columns in supernode J left to j:
  for (int k=supernodes[J];k<j;k++)
  {
    int jk = colpointers[k] + (j-k);
    int ik = jk;
    double& Fjk = F[jk];
    const double& Ljk = L[jk];
    for (int ij=s; ij<e; ij++)
    {
      // F[ik] = F[ik] - F[ij]*L[jk];
      // F[jk] = F[jk] - F[ij]*L[ik];
      F[ik] -= F[ij]*Ljk;
      Fjk   -= F[ij]*L[ik];
      ik++;
    }
  }
}

// Adjust column j for all columns in supernode K:
void ADcmod2(NumericVector& F,
            const NumericVector& L, int j, int K, int sz,
           NumericVector& t,
           const IntegerVector& indmap,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  // t is dense version of L[j], updated values at end of function:
  int i=0;
  for (int r = rowpointers[K+1] - 1;r>=rowpointers[K];r--)
  {
    int ndx = rowindices[r];
    int pos = colpointers[j+1] - 1 - indmap[ndx];
    t[i++] = F[pos];
  }

  // for all columns k in supernode K:
  for (int k=supernodes[K]; k<supernodes[K+1]; k++)
  {
    int jk = colpointers[k+1]-sz;
    int ik = jk;
    const double& Ljk = L[jk];
    double& Fjk = F[jk];
    for (int i=sz-1;i>=0;i--)
    {
      // F[ik] = F[ik] - F_ij*L[jk];
      // F[jk] = F[jk] - F_ij*L[ik];
      double F_ij = t[i];
      F[ik] -= F_ij*Ljk;
      Fjk   -= F_ij*L[ik];
      ik++;
    }
  }
}

void ADcdiv(NumericVector& F,
            const NumericVector& L, int j, const IntegerVector& colpointers)
{
  const int s = colpointers[j];
  const int e = colpointers[j+1];

  // update AD for column j:
  const double& Ls = L[s];
  double& Fs = F[s];
  for (int i = s + 1; i < e; i++)
  {
    // F[i] = F[i]/L[s];
    // F[s] = F[s] - L[i]*F[i];
    F[i] /= Ls;
    Fs -= L[i]*F[i];
  }
  //F[s] = Fs;
  F[s] = 0.5*F[s]/Ls;
}

void ADcholesky(NumericVector& F,
              const NumericVector& L,
              const IntegerVector& supernodes,
              const IntegerVector& rowpointers,
              const IntegerVector& colpointers,
              const IntegerVector& rowindices)
{
  const int N = colpointers.size() - 1;
  const int Nsupernodes = supernodes.size()-1;

  // linked lists, see section 4.2 Ng and Peyton
  IntegerVector HEAD(N,-1);
  IntegerVector LINK(Nsupernodes,-1);

  IntegerVector colhead = clone(rowpointers);
  IntegerVector coltop = clone(rowpointers);
  for (int J=0; J<Nsupernodes;J++)
  {
    int szNode = supernodes[J+1] - supernodes[J];
    coltop[J] += szNode-1;
    colhead[J] = rowpointers[J+1]-1;
    if (colhead[J] > coltop[J])
    {
      int rNdx = rowindices[colhead[J]];
      insert(HEAD, LINK, rNdx, J);
    }
  }
  IntegerVector indmap(N,0);
  NumericVector t(N);
  for (int J=Nsupernodes-1; J>=0;J--)
  {
    makeIndMap(indmap, J, rowpointers, rowindices);
    for (int j = supernodes[J+1]-1; j>=supernodes[J]; j--)
    {
      ADcdiv(F, L, j, colpointers);
      ADcmod1(F, L, j, J, supernodes, colpointers);

      int K = HEAD[j];
      while (K!=-1)
      {
        int nextK = LINK[K];
        colhead[K]--;
        if (colhead[K] > coltop[K])
        {
           int rNdx = rowindices[colhead[K]];
           insert(HEAD, LINK, rNdx, K);
        }
        int sz = rowpointers[K+1] - 1 - colhead[K];
        ADcmod2(F, L, j, K, sz, t, indmap, supernodes, rowpointers,colpointers,rowindices);
        K = nextK;
      }
      HEAD[j] = -1;
    }
  }
  return;
}

void initAD(NumericVector& F, const NumericVector& L, const IntegerVector& colpointers)
{
  std::fill(F.begin(), F.end(), 0.0);
  const int N = colpointers.size() - 1;
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    F[s] = 2.0/L[s];
  }
}

