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
  const double *l = L.begin();
  double *f = F.begin();

  int s = colpointers[j];
  int e = colpointers[j+1];
  // for all columns in supernode J left to j:
  for (int k=supernodes[J];k<j;k++)
  {
    int jk = colpointers[k] + (j-k);
    int ik = jk;
    double& fjk = f[jk];
    const double Ljk = l[jk];
    for (int ij=s; ij<e; ij++)
    {
      // F[ik] = F[ik] - F[ij]*L[jk];
      // F[jk] = F[jk] - F[ij]*L[ik];
      f[ik] -= f[ij]*Ljk;
      fjk   -= f[ij]*l[ik];
      ik++;
    }
  }
}


void ADcdiv(NumericVector& F,
            const NumericVector& L, int j, const IntegerVector& colpointers)
{
  const double *l = L.begin();
  double *f = F.begin();

  const int s = colpointers[j];
  const int e = colpointers[j+1];

  // update AD for column j:
  const double Ls = l[s];
  double Fs = f[s];
  for (int i = s + 1; i < e; i++)
  {
    // F[i] = F[i]/L[s];
    // F[s] = F[s] - L[i]*F[i];
    f[i] /= Ls;
    Fs -= l[i]*f[i];
  }
  //F[s] = Fs;
  f[s] = 0.5*Fs/Ls;
}

// ------------------------------------------------------------
// Reverse/AD update from source supernode K to target
// supernode J.
//
// done   = number of rows already consumed from the bottom
//          of K's off-diagonal row list.
//
// ncolup = number of rows currently belonging to J.
//
// The target rows are processed from bottom to top, as in the
// original scalar ADcholesky().
// ------------------------------------------------------------

void ADcmod2_sup(
    NumericVector& F,
    const NumericVector& L,
    int K,
    int done,
    int ncolup,
    NumericVector& t,
    const IntegerVector& indmap,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  if (ncolup <= 0)
    return;

  const double* l = L.begin();
  double* f = F.begin();
  double* tp = t.begin();

  const int row0 = rowpointers[K];
  const int eK   = rowpointers[K + 1];

  const int sCol = supernodes[K];
  const int eCol = supernodes[K + 1];

  // ----------------------------------------------------------
  // Process the affected target columns in reverse order.
  // ----------------------------------------------------------

  for (int p = 0; p < ncolup; ++p)
  {
    const int rj = eK - 1 - done - p;
    const int j = rowindices[rj];
    const int sz = done + p + 1;

    // --------------------------------------------------------
    // Gather exactly as in the original ADcmod2().
    // --------------------------------------------------------

    int i = 0;

    for (int r = eK - 1; r >= row0; --r)
    {
      const int ndx = rowindices[r];
      const int pos = colpointers[j + 1] - 1 - indmap[ndx];

      tp[i++] = f[pos];
    }

    // --------------------------------------------------------
    // Reverse update through all columns of K.
    // --------------------------------------------------------

    for (int k = sCol; k < eCol; ++k)
    {
      const int jk = colpointers[k + 1] - sz;
      int ik = jk;

      const double Ljk = l[jk];
      double& Fjk = f[jk];

      for (int i = sz - 1; i >= 0; --i)
      {
        const double Fij = tp[i];
        f[ik] -= Fij * Ljk;
        Fjk   -= Fij * l[ik];
        ++ik;
      }
    }
  }
}

void ADcholesky(
    NumericVector& F,
    const NumericVector& L,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  const int N = colpointers.size() - 1;
  const int Nsupernodes = supernodes.size() - 1;

  // SNODE[j] = supernode containing scalar row/column j
  IntegerVector SNODE(N);

  for (int J = 0; J < Nsupernodes; ++J)
  {
    for (int j = supernodes[J]; j < supernodes[J + 1];++j)
    {
      SNODE[j] = J;
    }
  }

  // ------------------------------------------------------------
  // Reverse supernodal linked lists.
  //
  // HEAD[J] = first source supernode K waiting to update J
  // LINK[K] = next source supernode on that list
  //
  // LENGTH[K] = number of rows of K already consumed from the
  //             bottom of its off-diagonal part.
  // ------------------------------------------------------------

  IntegerVector HEAD(Nsupernodes, -1);
  IntegerVector LINK(Nsupernodes, -1);
  IntegerVector LENGTH(Nsupernodes, 0);

  // ------------------------------------------------------------
  // Initially schedule every supernode on the list of the
  // supernode containing its last off-diagonal row.
  //
  // Unlike the forward case, HEAD and LINK are separate.
  // This is essential in the reverse traversal.
  // ------------------------------------------------------------

  for (int K = 0; K < Nsupernodes;++K)
  {
    const int width = supernodes[K + 1] - supernodes[K];
    const int len = rowpointers[K + 1] - rowpointers[K];
    const int offdiag = len - width;

    LENGTH[K] = 0;

    if (offdiag > 0)
    {
      const int last_row = rowindices[rowpointers[K + 1] - 1];
      const int J = SNODE[last_row];

      LINK[K] = HEAD[J];
      HEAD[J] = K;
    }
  }

  // Workspace
  IntegerVector indmap(N, 0);
  NumericVector t(N);

  // ------------------------------------------------------------
  // Reverse through supernodes
  // ------------------------------------------------------------

  for (int J = Nsupernodes - 1; J >= 0; --J)
  {
    const int j0 = supernodes[J];
    const int j1 = supernodes[J + 1];

    // Construct map for target supernode J.
    makeIndMap(indmap, J, rowpointers, rowindices);

    // Reverse operations internal to J.
    for (int j = j1 - 1;j >= j0;--j)
    {
      ADcdiv(F, L, j, colpointers);
      ADcmod1(F, L, j, J, supernodes,colpointers);
    }

    // Process all source supernodes K waiting to update J.
    int K = HEAD[J];
    HEAD[J] = -1;

    while (K != -1)
    {
      // Save next source supernode.
      const int nextK = LINK[K];
      const int done = LENGTH[K];
      const int row0 = rowpointers[K];
      const int eK =rowpointers[K + 1];
      const int widthK =supernodes[K + 1] - supernodes[K];
      const int lenK = eK - row0;
      const int offdiagK = lenK - widthK;

      // --------------------------------------------------------
      // Count the rows of K belonging to J, moving upward from
      // the current bottom.
      // --------------------------------------------------------

      int ncolup = 0;

      while (done + ncolup < offdiagK)
      {
        const int row = rowindices[eK - 1 - done - ncolup];

        if (row < j0)
          break;

        ++ncolup;
      }

      // Numerical reverse update.
      ADcmod2_sup(F, L, K, done, ncolup, t, indmap,
        supernodes, rowpointers, colpointers, rowindices);

      // Advance K further upward in its row list.
      const int newdone = done + ncolup;

      if (newdone < offdiagK)
      {
        const int next_row =rowindices[eK - 1 - newdone];
        const int nextJ = SNODE[next_row];
        LENGTH[K] = newdone;
        LINK[K] = HEAD[nextJ];
        HEAD[nextJ] = K;
      }
      else
      {
        LENGTH[K] = newdone;
        LINK[K] = -1;
      }
      K = nextK;
    }
  }
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

