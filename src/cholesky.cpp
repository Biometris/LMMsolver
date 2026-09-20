// sparse left-looking Cholesky using supernodes.
//
// For details see:
// Ng, Esmond G., and Barry W. Peyton.,
// "Block sparse Cholesky algorithms on advanced uniprocessor computers."
// SIAM Journal on Scientific Computing 14, no. 5 (1993): 1034-1056.
//
// Furrer, Reinhard, and Stephan R. Sain.
// "spam: A sparse matrix R package with emphasis on MCMC
// methods for Gaussian Markov random fields."
// Journal of Statistical Software 36 (2010): 1-25.

#include <Rcpp.h>
#include <set>
#include <vector>
#include "AuxFun.h"
#include "SparseMatrix.h"
#include "cholesky.h"
#include <chrono>

using namespace std::chrono;
using namespace Rcpp;
using namespace std;


// make indmap for supernode J:
void makeIndMap(IntegerVector& indmap,
                int J,
                const IntegerVector& rowpointers,
                const IntegerVector& rowindices)
{
  int s = rowpointers[J];
  int e = rowpointers[J+1];
  int l = 0;
  for (int i=e-1; i>=s;i--)
  {
    indmap[rowindices[i]] = l++;
  }
}

// Updates one target column j using 4 consecutive source
// columns k0,...,k0+3.
//
// For each row i:
//
//   y = L[i,j]
//   y -= L[i,k0  ] * L[j,k0  ]
//   y -= L[i,k0+1] * L[j,k0+1]
//   y -= L[i,k0+2] * L[j,k0+2]
//   y -= L[i,k0+3] * L[j,k0+3]
//   L[i,j] = y
//
// The target value is therefore loaded/stored only once.
// ============================================================

inline void update_column_cmod1_unroll4(
    double* l,
    int s,
    int e,
    int j,
    int k0,
    const IntegerVector& colpointers)
{
  // Starting positions of the active parts of the
  // four source columns.
  const double* p0 = l + colpointers[k0]     + (j - k0);
  const double* p1 = l + colpointers[k0 + 1] + (j - (k0 + 1));
  const double* p2 = l + colpointers[k0 + 2] + (j - (k0 + 2));
  const double* p3 = l + colpointers[k0 + 3] + (j - (k0 + 3));

  // The first element of each active source column is
  // L[j,k].
  const double a0 = *p0;
  const double a1 = *p1;
  const double a2 = *p2;
  const double a3 = *p3;

  // Update target column.
  for (int i = s; i < e; ++i)
  {
    double y = l[i];

    y -= (*p0++) * a0;
    y -= (*p1++) * a1;
    y -= (*p2++) * a2;
    y -= (*p3++) * a3;

    l[i] = y;
  }
}

// ============================================================
// cmod1 using source-column unrolling
// ============================================================

void cmod1(
    NumericVector& L,
    int j,
    int J,
    const IntegerVector& supernodes,
    const IntegerVector& colpointers)
{
  double* l = L.begin();

  const int s = colpointers[j];
  const int e = colpointers[j + 1];

  const int kstart = supernodes[J];
  const int nsrc = j - kstart;

  // Number of complete groups of 8.
  const int n4 =  (nsrc / 4) * 4;

  int k = kstart;

  // ----------------------------------------------------------
  // Groups of 4 source columns.
  // ----------------------------------------------------------

  const int kend = kstart + n4;

  for (; k < kend; k += 4)
  {
    update_column_cmod1_unroll4(l, s, e, j, k, colpointers);
  }

  // Remaining source columns, as in original cmod1
  for (; k < j; ++k)
  {
    const int jk = colpointers[k] + (j - k);

    int ik = jk;
    const double Ljk = l[jk];

    for (int i = s; i < e; ++i)
      l[i] -= l[ik++] * Ljk;
  }
}

void cdiv(NumericVector& L, int j, const IntegerVector& colpointers)
{
  double *l = L.begin();
  const int s = colpointers[j];
  const int e = colpointers[j+1];

  // pivot:
  l[s] = sqrt(l[s]);
  // update column j:
  double Ls = l[s];
  for (int i = s + 1; i < e; i++)
  {
    l[i] /= Ls;
  }
}


inline void update_column_cmod2_unroll4(
    double* l,
    double* t,
    int sz,
    int k0,
    const IntegerVector& colpointers)
{
  const double* p0 = l + colpointers[k0 + 1] - sz;
  const double* p1 = l + colpointers[k0 + 2] - sz;
  const double* p2 = l + colpointers[k0 + 3] - sz;
  const double* p3 = l + colpointers[k0 + 4] - sz;

  const double a0 = *p0;
  const double a1 = *p1;
  const double a2 = *p2;
  const double a3 = *p3;

  double* tptr = t + sz - 1;

  for (int i = 0; i < sz; ++i)
  {
    *tptr-- += (*p0++) * a0
             + (*p1++) * a1
             + (*p2++) * a2
             + (*p3++) * a3;
  }
}

// cmod2 using 4 source-column unrolling
void cmod2(
    NumericVector& L,
    int J,
    int K,
    int khead,
    int klen,
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

  double* l = L.begin();
  double* tp = t.begin();

  const int eK = rowpointers[K + 1];
  const int sCol = supernodes[K];
  const int eCol = supernodes[K + 1];
  const int srcWidth = eCol - sCol;

  // Four source columns at a time.
  const int n4 = (srcWidth / 4) * 4;

  // One target column at a time.
  for (int p = 0; p < ncolup; ++p)
  {
    const int j = rowindices[khead + p];
    const int sz = klen - p;

    // Special case: one source column.
    if (srcWidth == 1)
    {
      const int src_end = colpointers[sCol + 1];
      const int jk = src_end - sz;
      const double Ljk = l[jk];

      const double* lptr = l + src_end - 1;

      int r = eK - 1;
      const int ref_pos = colpointers[j + 1] - 1;

      for (int i = 0; i < sz; ++i)
      {
        const int ndx = rowindices[r--];
        const int pos = ref_pos - indmap[ndx];

        l[pos] -= *lptr-- * Ljk;
      }

      continue;
    }

    // Initialise t.
    for (int i = 0; i < sz; ++i)
      tp[i] = 0.0;

    int k = sCol;

    for (; k < sCol + n4; k += 4)
    {
      update_column_cmod2_unroll4(l, tp, sz, k, colpointers);
    }

    // Remaining source columns.
    for (; k < eCol; ++k)
    {
      const int jk = colpointers[k + 1] - sz;
      const double Ljk = l[jk];
      const double* lptr = l + jk;
      double* tptr = tp + sz - 1;

      for (int i = 0; i < sz; ++i)
      {
        *tptr-- += *lptr++ * Ljk;
      }
    }

    // Scatter.
    int r = eK - 1;

    const int ref_pos = colpointers[j + 1] - 1;

    for (int i = 0; i < sz; ++i)
    {
      const int ndx = rowindices[r--];
      const int pos = ref_pos - indmap[ndx];

      l[pos] -= tp[i];
    }
  }
}


void cholesky(
    NumericVector& L,
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
  // Supernodal linked lists
  //
  // LINK[J]   = head of list of source supernodes waiting
  //             to update J
  //
  // LENGTH[K] = active suffix length of source supernode K
  // ------------------------------------------------------------

  IntegerVector LINK(Nsupernodes, -1);
  IntegerVector LENGTH(Nsupernodes, 0);

  // ------------------------------------------------------------
  // Workspace
  // ------------------------------------------------------------

  IntegerVector indmap(N, 0);

  // Temporary vector used by cmod2.
  NumericVector t(N);

  // ------------------------------------------------------------
  // Process supernodes in order
  // ------------------------------------------------------------

  for (int J = 0; J < Nsupernodes; ++J)
  {
    const int j0 = supernodes[J];
    const int j1 = supernodes[J + 1];

    makeIndMap(indmap, J, rowpointers, rowindices);

    // Process all source supernodes currently waiting for J.
    int K = LINK[J];

    LINK[J] = -1;

    while (K != -1)
    {
      const int nextK = LINK[K];
      const int klen = LENGTH[K];

      // First active row of K.
      const int khead = rowpointers[K + 1] - klen;

      // --------------------------------------------------------
      // Determine how many active rows of K belong to J.
      // --------------------------------------------------------

      int ncolup = 0;

      while (ncolup < klen &&
             rowindices[khead + ncolup] < j1)
      {
        ++ncolup;
      }

      // --------------------------------------------------------
      // Numerical update.
      //
      // One target column at a time, with source columns of K
      // processed four at a time.
      // --------------------------------------------------------

      if (ncolup > 0)
      {
        cmod2(L, J, K, khead, klen, ncolup, t, indmap,
          supernodes, rowpointers, colpointers, rowindices);
      }

      // --------------------------------------------------------
      // K may still have an active suffix.
      //
      // If so, put K on the list of the next supernode it
      // contributes to.
      // --------------------------------------------------------

      if (klen > ncolup)
      {
        const int next_row = rowindices[khead + ncolup];
        const int nextJ = SNODE[next_row];

        LENGTH[K] = klen - ncolup;

        LINK[K] = LINK[nextJ];
        LINK[nextJ] = K;
      }
      else
      {
        LENGTH[K] = 0;
        LINK[K] = -1;
      }

      K =
        nextK;
    }

    // ----------------------------------------------------------
    // Factor supernode J
    // ----------------------------------------------------------

    for (int j = j0; j < j1; ++j)
    {
      cmod1(L, j, J, supernodes, colpointers);
      cdiv(L, j, colpointers);
    }

    // ----------------------------------------------------------
    // Schedule J's own update.
    //
    // The first 'width' entries of J's row list correspond to
    // its diagonal supernode block. Everything after that is
    // the update that J still has to deliver.
    // ----------------------------------------------------------

    const int width = j1 - j0;

    const int len = rowpointers[J + 1] - rowpointers[J];

    LENGTH[J] = len - width;

    if (LENGTH[J] > 0)
    {
      const int next_row = rowindices[rowpointers[J] + width];
      const int nextJ = SNODE[next_row];

      LINK[J] = LINK[nextJ];
      LINK[nextJ] = J;
    }
    else
    {
      LENGTH[J] = 0;
      LINK[J] = -1;
    }
  }
}


double logdet(const NumericVector& L, const IntegerVector& colpointers)
{
  const int N = colpointers.size() - 1;
  double sum = 0;
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    sum += 2.0*log(L[s]);
  }
  return sum;
}


NumericVector forwardCholesky(
    const NumericVector& L,
    const NumericVector& b,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices,
    const IntegerVector& pivot,
    const IntegerVector& invpivot)
{
  const int Nsupernodes = supernodes.size()-1;
  const int N = colpointers.size() - 1;
  NumericVector x(N);
  NumericVector Pb(N);  // permutation of P
  NumericVector sum(N); // sum for each column
  for (int i=0;i<N;i++)
  {
    Pb[i] = b[pivot[i]];
  }
  for (int J=0; J<Nsupernodes;J++)
  {
    int s = rowpointers[J];
    for (int j=supernodes[J]; j<supernodes[J+1]; j++)
    {
      const double x_j = (Pb[j]-sum[j])/L[colpointers[j]];
      x[j] = x_j;
      // the non-diagonal elements of column j
      for (int ndx = colpointers[j]+1, k=s+1; ndx < colpointers[j+1]; ndx++)
      {
        int i = rowindices[k++];
        sum[i] += L[ndx]*x_j;
      }
      s++;
    }
  }

  NumericVector xP(N); // inverse permutation
  for (int i=0;i<N;i++) {
    xP[i] = x[invpivot[i]];
  }
  return xP;

}


NumericVector backwardCholesky(
    const NumericVector& L,
    const NumericVector& b,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices,
    const IntegerVector& pivot,
    const IntegerVector& invpivot)
{
  const int Nsupernodes = supernodes.size()-1;
  const int N = colpointers.size() - 1;
  NumericVector x(N);
  NumericVector Pb(N);  // permutation of P
  NumericVector sum(N); // sum for each column
  for (int i=0;i<N;i++)
  {
    Pb[i] = b[pivot[i]];
  }
  for (int J=Nsupernodes-1; J>=0;J--)
  {
    int NodeSz = supernodes[J+1] - supernodes[J];
    int s = rowpointers[J]+(NodeSz-1);
    for (int j=supernodes[J+1]-1; j>=supernodes[J]; j--)
    {
      //Rcout << "column " << j << endl;
      double alpha = L[colpointers[j]];
      double x_j = Pb[j];
      // the non-diagonal elements of column j
      for (int ndx = colpointers[j]+1, k=s+1; ndx < colpointers[j+1]; ndx++)
      {
        int i = rowindices[k++];
        x_j -= L[ndx]*x[i];
      }
      x[j] = x_j/alpha;
      s--;
    }
  }

  NumericVector xP(N); // inverse permutation
  for (int i=0;i<N;i++) {
    xP[i] = x[invpivot[i]];
  }
  return xP;
}

