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

// j is current column in Supernode J
void cmod1(NumericVector& L, int j, int J,
           const IntegerVector& supernodes,
           const IntegerVector& colpointers)
{
  double *l = L.begin();
  const int s = colpointers[j];
  const int e = colpointers[j+1];
  // for all columns in supernode J left to j:
  for (int k=supernodes[J];k<j;k++)
  {
    const int jk = colpointers[k] + (j-k);
    int ik = jk;
    const double Ljk = l[jk];
    for (int ij=s; ij<e; ij++)
    {
       l[ij] -= l[ik++]*Ljk;
       //ik++;
    }
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

// ============================================================
// Factor and solve a supernode of k columns.
//
// The current column block contains:
//
//     [ A11 ]
//     [ A21 ]
//
// where A11 is k x k symmetric positive definite.
//
// Computes:
//
//     A11 = L11 L11^T
//
// followed by:
//
//     L21 = A21 L11^{-T}
//
// This replaces the sequence of cmod1() + cdiv() for the
// columns within one supernode.
//
// No special optimisation; intended for small k (e.g. 4 or 8).
// ============================================================

void cdiv_block(
    NumericVector& L,
    int j0,
    int k,
    const IntegerVector& colpointers)
{
  double* l = L.begin();

  // ----------------------------------------------------------
  // 1. Cholesky factorisation of the k x k diagonal block.
  //
  // Column c starts at colpointers[j0 + c].
  // Within that column, row r (r >= c) is at offset r-c.
  // ----------------------------------------------------------

  for (int c = 0; c < k; ++c)
  {
    const int sc =
      colpointers[j0 + c];

    // --------------------------------------------------------
    // Diagonal element.
    // --------------------------------------------------------

    double d =
      l[sc];

    for (int h = 0; h < c; ++h)
    {
      const int sh =
        colpointers[j0 + h];

      const double Lch =
        l[sh + (c - h)];

      d -=
        Lch * Lch;
    }

    l[sc] =
      sqrt(d);

    const double Lcc =
      l[sc];

    // --------------------------------------------------------
    // Elements below the diagonal within the supernode.
    // --------------------------------------------------------

    for (int r = c + 1; r < k; ++r)
    {
      const int pos =
        sc + (r - c);

      double x =
        l[pos];

      for (int h = 0; h < c; ++h)
      {
        const int sh =
          colpointers[j0 + h];

        const double Lrh =
          l[sh + (r - h)];

        const double Lch =
          l[sh + (c - h)];

        x -=
          Lrh * Lch;
      }

      l[pos] =
        x / Lcc;
    }
  }


  // ----------------------------------------------------------
  // 2. Forward/triangular solve for the off-diagonal block.
  //
  // The first k entries of each column belong to A11.
  // The remaining entries form A21.
  //
  // The off-diagonal row pattern is common to the columns
  // of the supernode.
  // ----------------------------------------------------------

  const int s0 =
    colpointers[j0];

  const int e0 =
    colpointers[j0 + 1];

  const int noffdiag =
    e0 - s0 - k;

  for (int t = 0; t < noffdiag; ++t)
  {
    for (int c = 0; c < k; ++c)
    {
      const int jc =
        j0 + c;

      const int pos =
        colpointers[jc] + (k - c) + t;

      double x =
        l[pos];

      for (int h = 0; h < c; ++h)
      {
        const int jh =
          j0 + h;

        const double Lrh =
          l[colpointers[jh] + (k - h) + t];

        const double Lch =
          l[colpointers[jh] + (c - h)];

        x -=
          Lrh * Lch;
      }

      const double Lcc =
        l[colpointers[jc]];

      l[pos] =
        x / Lcc;
    }
  }
}

// ============================================================
// Standard cmod2 update
//
// Updates ncolup target columns of J from source supernode K
// using the original t-based implementation.
// ============================================================

inline void cmod2_default(
    double* l,
    double* tp,
    int K,
    int khead,
    int klen,
    int ncolup,
    int eK,
    int sCol,
    int srcWidth,
    const IntegerVector& indmap,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  for (int p = 0; p < ncolup; ++p)
  {
    const int j  = rowindices[khead + p];
    const int sz = klen - p;

    // --------------------------------------------------------
    // Initialise t
    // --------------------------------------------------------

    for (int i = 0; i < sz; ++i)
      tp[i] = 0.0;

    // --------------------------------------------------------
    // Accumulate source supernode contribution
    // --------------------------------------------------------

    for (int k = 0; k < srcWidth; ++k)
    {
      const int jk =
        colpointers[sCol + k + 1] - sz;

      const double Ljk =
        l[jk];

      double* tptr =
        tp + sz - 1;

      double* lptr =
        l + jk;

      for (int i = 0; i < sz; ++i)
        *tptr-- += *lptr++ * Ljk;
    }

    // --------------------------------------------------------
    // Scatter
    // --------------------------------------------------------

    int r = eK - 1;

    const int ref_pos =
      colpointers[j + 1] - 1;

    double* tptr = tp;

    for (int i = 0; i < sz; ++i)
    {
      const int ndx =
        rowindices[r--];

      const int pos =
        ref_pos - indmap[ndx];

      l[pos] -= *tptr++;
    }
  }
}

// ============================================================
// cmod2_sup
// ============================================================

void cmod2_sup(
    NumericVector& L,
    int J,
    int K,
    int khead,
    int klen,
    int ncolup,
    NumericVector& t,
    std::vector<double>& A,
    std::vector<double>& Bt,
    const IntegerVector& indmap,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  if (ncolup <= 0)
    return;

  double* l  = L.begin();
  double* tp = t.begin();

  const int eK =
    rowpointers[K + 1];

  const int sCol =
    supernodes[K];

  const int eCol =
    supernodes[K + 1];

  const int srcWidth =
    eCol - sCol;

  // ==========================================================
  // Small source supernodes:
  // use standard cmod2.
  // ==========================================================

  if (srcWidth <= 3)
  {
    cmod2_default(
      l,
      tp,
      K,
      khead,
      klen,
      ncolup,
      eK,
      sCol,
      srcWidth,
      indmap,
      colpointers,
      rowindices
    );

    return;
  }


  // ==========================================================
  // Source-column base positions
  // ==========================================================

  std::vector<int> srcBase(srcWidth);

  for (int k = 0; k < srcWidth; ++k)
  {
    srcBase[k] =
      colpointers[sCol + k + 1] - klen;
  }


  // ==========================================================
  // Process target columns in groups of four.
  // ==========================================================

  const int ncol4 =
    (ncolup / 4) * 4;

  int p0 = 0;

  for (; p0 < ncol4; p0 += 4)
  {
    // --------------------------------------------------------
    // Check whether the target block overlaps source K.
    // If so, use the original in-place update.
    // --------------------------------------------------------

    bool useOptimized = true;

    for (int j = 0; j < 4; ++j)
    {
      const int target =
        rowindices[khead + p0 + j];

      if (target >= sCol && target < eCol)
      {
        useOptimized = false;
        break;
      }
    }

    if (!useOptimized)
    {
      cmod2_default(
        l,
        tp,
        K,
        khead + p0,
        klen,
        4,
        eK,
        sCol,
        srcWidth,
        indmap,
        colpointers,
        rowindices
      );

      continue;
    }


    // ========================================================
    // Build Bt for this target block
    // ========================================================

    for (int k = 0; k < srcWidth; ++k)
    {
      const int base =
        srcBase[k];

      Bt[k * 4 + 0] =
        l[base + p0 + 0];

      Bt[k * 4 + 1] =
        l[base + p0 + 1];

      Bt[k * 4 + 2] =
        l[base + p0 + 2];

      Bt[k * 4 + 3] =
        l[base + p0 + 3];
    }


    // ========================================================
    // Row blocks
    // ========================================================

    int r0 = p0;

    for (; r0 + 4 <= klen; r0 += 4)
    {
      update_block<4,4>(
          l,
          A,
          Bt,
          srcWidth,
          r0,
          p0,
          klen,
          eK,
          srcBase,
          khead,
          indmap,
          colpointers,
          rowindices
      );
    }


    // --------------------------------------------------------
    // Tail rows
    // --------------------------------------------------------

    if (r0 < klen)
    {
      const int M =
        klen - r0;

      if (M == 1)
      {
        update_block<1,4>(
            l,
            A,
            Bt,
            srcWidth,
            r0,
            p0,
            klen,
            eK,
            srcBase,
            khead,
            indmap,
            colpointers,
            rowindices
        );
      }
      else if (M == 2)
      {
        update_block<2,4>(
            l,
            A,
            Bt,
            srcWidth,
            r0,
            p0,
            klen,
            eK,
            srcBase,
            khead,
            indmap,
            colpointers,
            rowindices
        );
      }
      else
      {
        update_block<3,4>(
            l,
            A,
            Bt,
            srcWidth,
            r0,
            p0,
            klen,
            eK,
            srcBase,
            khead,
            indmap,
            colpointers,
            rowindices
        );
      }
    }
  }


  // ==========================================================
  // Remaining target columns
  // ==========================================================

  if (p0 < ncolup)
  {
    cmod2_default(
      l,
      tp,
      K,
      khead + p0,
      klen - p0,
      ncolup - p0,
      eK,
      sCol,
      srcWidth,
      indmap,
      colpointers,
      rowindices
    );
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

  // SNODE[j] = supernode containing scalar row/column j
  IntegerVector SNODE(N);

  for (int J = 0; J < Nsupernodes; ++J)
  {
    for (int j = supernodes[J]; j < supernodes[J + 1]; ++j)
      SNODE[j] = J;
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

  int maxSrcWidth = 0;

  for (int K = 0; K < supernodes.size() - 1; ++K)
    maxSrcWidth = std::max(
      maxSrcWidth,
      supernodes[K + 1] - supernodes[K]
    );

  NumericVector t(N);
  std::vector<double> A(4 * maxSrcWidth);
  std::vector<double> Bt(4 * maxSrcWidth);

  // ------------------------------------------------------------
  // Process supernodes in order
  // ------------------------------------------------------------

  for (int J = 0; J < Nsupernodes; ++J)
  {
    const int j0 = supernodes[J];
    const int j1 = supernodes[J + 1];

    makeIndMap(indmap, J, rowpointers, rowindices);

    // ----------------------------------------------------------
    // Process all source supernodes currently waiting for J.
    //
    // Important: LINK[K] is the next pointer when K is on
    // another list, so save it before changing LINK[K].
    // ----------------------------------------------------------

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
      //
      // Because K was put on J's list when its first active row
      // was in J, rowindices[khead] should be >= j0.
      // --------------------------------------------------------

      int ncolup = 0;

      while (ncolup < klen && rowindices[khead + ncolup] < j1)
      {
        ++ncolup;
      }

      // Numerical update.
      cmod2_sup(L, J, K, khead, klen, ncolup, t, A, Bt, indmap,
        supernodes, rowpointers, colpointers, rowindices
      );

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

      K = nextK;
    }

    // ----------------------------------------------------------
    // Phase 2:
    // factor supernode J
    // ----------------------------------------------------------

    for (int j = j0; j < j1; ++j)
    {
      cmod1(L,j,J, supernodes, colpointers);
      cdiv_block(L, j, 1, colpointers);
    }

    // ----------------------------------------------------------
    // Now schedule J's own update.
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
      // First row below the diagonal block.
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

