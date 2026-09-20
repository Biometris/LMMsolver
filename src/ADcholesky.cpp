// Backwards Automatic Differentiation of the Cholesky Algorithm
// for calculating partial derivatives of the log-determinant of
// positive definite symmetric sparse matrices.
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
//
// Implementation by Martin Boer, 2026.

#include <Rcpp.h>
#include <set>
#include <vector>
#include "AuxFun.h"
#include "SparseMatrix.h"
#include "cholesky.h"
#include "ADcholesky.h"

using namespace Rcpp;
using namespace std;


// ============================================================
// AD update of one target column using 4 consecutive source
// columns.
//
// For each source column k:
//
//   f[i,k] -= f[i,j] * L[j,k]
//   f[j,k] -= f[i,j] * L[i,k]
//
// The first source entry L[j,k] is also F[j,k], so the first
// row contributes twice to F[j,k]. This is handled explicitly.
// ============================================================

inline void update_AD_column_cmod1_unroll4(
    double* f,
    const double* l,
    int s,
    int e,
    int j,
    int k0,
    const IntegerVector& colpointers)
{
  // ----------------------------------------------------------
  // First positions of the four source columns.
  // These are the L[j,k] values and also the F[j,k] values.
  // ----------------------------------------------------------

  const int jk0 = colpointers[k0] + (j - k0);
  const int jk1 = colpointers[k0 + 1] + (j - (k0 + 1));
  const int jk2 = colpointers[k0 + 2] + (j - (k0 + 2));
  const int jk3 = colpointers[k0 + 3] + (j - (k0 + 3));

  const double a0 = l[jk0];
  const double a1 = l[jk1];
  const double a2 = l[jk2];
  const double a3 = l[jk3];

  // ----------------------------------------------------------
  // First target row.
  //
  // Here source F[j,k] and fjk are the same element, so the
  // original code subtracts the contribution twice.
  // ----------------------------------------------------------

  const double f0 = f[s];

  double g0 = f[jk0] - 2.0 * f0 * a0;
  double g1 = f[jk1] - 2.0 * f0 * a1;
  double g2 = f[jk2] - 2.0 * f0 * a2;
  double g3 = f[jk3] - 2.0 * f0 * a3;

  // Remaining rows.

  const double* p0 = l + jk0 + 1;
  const double* p1 = l + jk1 + 1;
  const double* p2 = l + jk2 + 1;
  const double* p3 = l + jk3 + 1;

  double* q0 = f + jk0 + 1;
  double* q1 = f + jk1 + 1;
  double* q2 = f + jk2 + 1;
  double* q3 = f + jk3 + 1;

  for (int i = s + 1; i < e; ++i)
  {
    const double fij = f[i];

    *q0++ -= fij * a0;
    g0   -= fij * (*p0++);

    *q1++ -= fij * a1;
    g1   -= fij * (*p1++);

    *q2++ -= fij * a2;
    g2   -= fij * (*p2++);

    *q3++ -= fij * a3;
    g3   -= fij * (*p3++);
  }

  // ----------------------------------------------------------
  // Write back F[j,k].
  // ----------------------------------------------------------

  f[jk0] = g0;
  f[jk1] = g1;
  f[jk2] = g2;
  f[jk3] = g3;
}


// ============================================================
// ADcmod1 using source-column unrolling
// ============================================================

void ADcmod1(
    NumericVector& F,
    const NumericVector& L,
    int j,
    int J,
    const IntegerVector& supernodes,
    const IntegerVector& colpointers)
{
  const double* l = L.begin();

  double* f = F.begin();

  const int s = colpointers[j];
  const int e = colpointers[j + 1];

  const int kstart = supernodes[J];
  const int nsrc = j - kstart;

  // Complete groups of four source columns.
  const int n4 = (nsrc / 4) * 4;

  int k = kstart;

  for (; k < kstart + n4; k += 4)
  {
    update_AD_column_cmod1_unroll4(f, l, s, e, j, k, colpointers);
  }

  // ----------------------------------------------------------
  // Remaining source columns.
  for (; k < j; ++k)
  {
    const int jk = colpointers[k] + (j - k);
    int ik = jk;

    double& fjk = f[jk];

    const double Ljk = l[jk];

    for (int ij = s; ij < e; ++ij)
    {
      f[ik] -= f[ij] * Ljk;
      fjk   -= f[ij] * l[ik];
      ++ik;
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

// ============================================================
// AD update of one target column from 4 source columns.
// ============================================================

inline void update_AD_column_cmod2_unroll4(
    double* f,
    const double* l,
    const double* t,
    int sz,
    int k0,
    const IntegerVector& colpointers)
{
  const int p0 = colpointers[k0 + 1] - sz;
  const int p1 = colpointers[k0 + 2] - sz;
  const int p2 = colpointers[k0 + 3] - sz;
  const int p3 = colpointers[k0 + 4] - sz;

  const double a0 = l[p0];
  const double a1 = l[p1];
  const double a2 = l[p2];
  const double a3 = l[p3];

  // ----------------------------------------------------------
  // F[j,k] = F[j,k]
  //             - 2 * t[sz-1] * L[j,k]
  //
  // This is the overlap of the two reverse updates.
  // ----------------------------------------------------------

  const double fij = t[sz - 1];

  double g0 = f[p0] - 2.0 * fij * a0;
  double g1 = f[p1] - 2.0 * fij * a1;
  double g2 = f[p2] - 2.0 * fij * a2;
  double g3 = f[p3] - 2.0 * fij * a3;

  // Remaining source elements.
  const double* lp0 = l + p0 + 1;
  const double* lp1 = l + p1 + 1;
  const double* lp2 = l + p2 + 1;
  const double* lp3 = l + p3 + 1;
  double* fp0 = f + p0 + 1;
  double* fp1 = f + p1 + 1;
  double* fp2 = f + p2 + 1;
  double* fp3 = f + p3 + 1;

  for (int i = sz - 2; i >= 0; --i)
  {
    const double Fij = t[i];

    *fp0 -= Fij * a0;
    *fp1 -= Fij * a1;
    *fp2 -= Fij * a2;
    *fp3 -= Fij * a3;

    g0 -= Fij * (*lp0);
    g1 -= Fij * (*lp1);
    g2 -= Fij * (*lp2);
    g3 -= Fij * (*lp3);

    ++fp0;
    ++fp1;
    ++fp2;
    ++fp3;

    ++lp0;
    ++lp1;
    ++lp2;
    ++lp3;
  }

  f[p0] = g0;
  f[p1] = g1;
  f[p2] = g2;
  f[p3] = g3;
}


void ADcmod2(
    NumericVector& F,
    const NumericVector& L,
    int j,
    int K,
    int sz,
    NumericVector& t,
    const IntegerVector& indmap,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  const double* l = L.begin();
  double* f = F.begin();
  double* tp = t.begin();

  const int eK = rowpointers[K + 1];
  const int sCol = supernodes[K];
  const int eCol = supernodes[K + 1];
  const int srcWidth = eCol - sCol;

  // Special case: one source column.
  if (srcWidth == 1)
  {
    const int jk = colpointers[sCol + 1] - sz;
    const double Ljk = l[jk];
    double& Fjk = f[jk];

    int r = eK - sz;
    int ik = jk;

    const int ref_pos = colpointers[j + 1] - 1;

    for (int i = sz - 1; i >= 0; --i)
    {
      const int ndx = rowindices[r++];

      const int pos = ref_pos - indmap[ndx];

      // Copy before updating f[ik], since ik may equal pos.
      const double Fij = f[pos];

      f[ik] -= Fij * Ljk;
      Fjk    -= Fij * l[ik];

      ++ik;
    }

    return;
  }

  // Gather target AD values.
  int r = eK - 1;

  const int ref_pos = colpointers[j + 1] - 1;

  for (int i = 0; i < sz; ++i)
  {
    const int ndx = rowindices[r--];

    const int pos = ref_pos - indmap[ndx];

    tp[i] = f[pos];
  }

  // Four source columns at a time.
  const int n4 = (srcWidth / 4) * 4;

  int k = sCol;

  for (; k < sCol + n4; k += 4)
  {
    update_AD_column_cmod2_unroll4(f, l, tp, sz, k, colpointers);
  }

  // Remaining source columns.
  for (; k < eCol; ++k)
  {
    const int jk = colpointers[k + 1] - sz;
    const double Ljk = l[jk];
    double& Fjk = f[jk];

    int ik = jk;

    for (int i = sz - 1; i >= 0; --i)
    {
      const double Fij = tp[i];

      f[ik] -= Fij * Ljk;
      Fjk    -= Fij * l[ik];

      ++ik;
    }
  }
}


// ADcholesky
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

  // ------------------------------------------------------------
  // SNODE[j] = supernode containing scalar row/column j
  // ------------------------------------------------------------

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
  // HEAD[J]   = first source supernode K waiting to update J
  // LINK[K]   = next source supernode on that list
  //
  // LENGTH[K] = number of rows of K already consumed from the
  //              bottom of its off-diagonal part.
  // ------------------------------------------------------------

  IntegerVector HEAD(Nsupernodes, -1);
  IntegerVector LINK(Nsupernodes, -1);
  IntegerVector LENGTH(Nsupernodes, 0);

  // ------------------------------------------------------------
  // Initially schedule every supernode on the list of the
  // supernode containing its last off-diagonal row.
  // ------------------------------------------------------------

  for (int K = 0; K < Nsupernodes; ++K)
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

  // ------------------------------------------------------------
  // Workspace
  // ------------------------------------------------------------

  IntegerVector indmap(N, 0);

  NumericVector t(N);

  // ------------------------------------------------------------
  // Reverse through supernodes
  // ------------------------------------------------------------

  for (int J = Nsupernodes - 1;J >= 0;--J)
  {
    const int j0 = supernodes[J];
    const int j1 = supernodes[J + 1];

    // Construct map for target supernode J.
    makeIndMap(indmap, J, rowpointers, rowindices);

    // Reverse operations internal to J.
    for (int j = j1 - 1; j >= j0;--j)
    {
      ADcdiv(F, L, j, colpointers);
      ADcmod1(F, L, j, J, supernodes, colpointers);
    }

    // ----------------------------------------------------------
    // Process all source supernodes K waiting to update J.
    // ----------------------------------------------------------

    int K = HEAD[J];

    HEAD[J] = -1;

    while (K != -1)
    {
      // --------------------------------------------------------
      // Save next source supernode.
      // --------------------------------------------------------

      const int nextK = LINK[K];
      const int done = LENGTH[K];

      const int row0 = rowpointers[K];
      const int eK = rowpointers[K + 1];

      const int widthK = supernodes[K + 1] - supernodes[K];

      const int lenK = eK - row0;

      const int offdiagK = lenK - widthK;

      // --------------------------------------------------------
      // Count rows of K belonging to J, moving upward from
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

      // --------------------------------------------------------
      // Reverse update, one target column at a time.
      for (int p = 0; p < ncolup;++p)
      {
        const int j = rowindices[eK - 1 - done - p];
        const int sz = done + p + 1;

        ADcmod2(F, L, j, K, sz, t, indmap, supernodes,
            rowpointers, colpointers, rowindices);
      }

      // --------------------------------------------------------
      // Advance K further upward in its row list.
      // --------------------------------------------------------

      const int newdone = done + ncolup;

      if (newdone < offdiagK)
      {
        const int next_row = rowindices[eK - 1 - newdone];

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

