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

/*
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
*/

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

inline void update_column_unroll4(
    double* l,
    int s,
    int e,
    int j,
    int k0,
    const IntegerVector& colpointers)
{
  // Starting positions of the active parts of the
  // eight source columns.
  double* p0 =
    l + colpointers[k0]     + (j - k0);

  double* p1 =
    l + colpointers[k0 + 1] + (j - (k0 + 1));

  double* p2 =
    l + colpointers[k0 + 2] + (j - (k0 + 2));

  double* p3 =
    l + colpointers[k0 + 3] + (j - (k0 + 3));


  // The first element of each active source column is
  // L[j,k].
  const double a0 = *p0;
  const double a1 = *p1;
  const double a2 = *p2;
  const double a3 = *p3;

  // ----------------------------------------------------------
  // Update target column.
  // ----------------------------------------------------------

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
  double* l =
    L.begin();

  const int s =
    colpointers[j];

  const int e =
    colpointers[j + 1];

  const int kstart =
    supernodes[J];

  const int nsrc =
    j - kstart;

  // Number of complete groups of 8.
  const int n4 =
    (nsrc / 4) * 4;

  int k =
    kstart;

  // ----------------------------------------------------------
  // Groups of 4 source columns.
  // ----------------------------------------------------------

  const int kend =
    kstart + n4;

  for (; k < kend; k += 4)
  {
    update_column_unroll4(
      l,
      s,
      e,
      j,
      k,
      colpointers
    );
  }

  // ----------------------------------------------------------
  // Remaining source columns.
  //
  // Keep original cmod1 code for the tail.
  // ----------------------------------------------------------

  for (; k < j; ++k)
  {
    const int jk =
      colpointers[k] + (j - k);

    int ik =
      jk;

    const double Ljk =
      l[jk];

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


void cholesky_wh_timer(
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
      cdiv(L, j, colpointers);
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


void cholesky(
    NumericVector& L,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  using clock = std::chrono::steady_clock;

  const int N =
    colpointers.size() - 1;

  const int Nsupernodes =
    supernodes.size() - 1;

  // ------------------------------------------------------------
  // Timing
  // ------------------------------------------------------------

  const auto tTotal0 = clock::now();

  double tIndMap = 0.0;
  double tCmod2  = 0.0;
  double tCmod1  = 0.0;
  double tCdiv   = 0.0;

  long long nCmod2 = 0;
  long long nCmod1 = 0;
  long long nCdiv  = 0;

  // ------------------------------------------------------------
  // SNODE[j] = supernode containing scalar row/column j
  // ------------------------------------------------------------

  IntegerVector SNODE(N);

  for (int J = 0; J < Nsupernodes; ++J)
  {
    for (int j = supernodes[J];
         j < supernodes[J + 1];
         ++j)
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

  int maxSrcWidth = 0;

  for (int K = 0; K < supernodes.size() - 1; ++K)
  {
    maxSrcWidth = std::max(
      maxSrcWidth,
      supernodes[K + 1] - supernodes[K]
    );
  }

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

    // ----------------------------------------------------------
    // Build index map
    // ----------------------------------------------------------

    {
      const auto t0 = clock::now();

      makeIndMap(
        indmap,
        J,
        rowpointers,
        rowindices
      );

      const auto t1 = clock::now();

      tIndMap +=
        std::chrono::duration<double>(t1 - t0).count();
    }

    // ----------------------------------------------------------
    // Process all source supernodes currently waiting for J.
    // ----------------------------------------------------------

    int K = LINK[J];
    LINK[J] = -1;

    while (K != -1)
    {
      const int nextK = LINK[K];
      const int klen = LENGTH[K];

      // First active row of K.
      const int khead =
        rowpointers[K + 1] - klen;

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
      // --------------------------------------------------------

      if (ncolup > 0)
      {
        const auto t0 = clock::now();

        cmod2_sup(
          L,
          J,
          K,
          khead,
          klen,
          ncolup,
          t,
          A,
          Bt,
          indmap,
          supernodes,
          rowpointers,
          colpointers,
          rowindices
        );

        const auto t1 = clock::now();

        tCmod2 +=
          std::chrono::duration<double>(t1 - t0).count();

        ++nCmod2;
      }

      // --------------------------------------------------------
      // K may still have an active suffix.
      // --------------------------------------------------------

      if (klen > ncolup)
      {
        const int next_row =
          rowindices[khead + ncolup];

        const int nextJ =
          SNODE[next_row];

        LENGTH[K] =
          klen - ncolup;

        LINK[K] =
          LINK[nextJ];

        LINK[nextJ] =
          K;
      }
      else
      {
        LENGTH[K] = 0;
        LINK[K] = -1;
      }

      K = nextK;
    }

    // ----------------------------------------------------------
    // Factor supernode J
    // ----------------------------------------------------------

    for (int j = j0; j < j1; ++j)
    {
      // --------------------------------------------------------
      // cmod1
      // --------------------------------------------------------

      {
        const auto t0 = clock::now();

        cmod1(
          L,
          j,
          J,
          supernodes,
          colpointers
        );

        const auto t1 = clock::now();

        tCmod1 +=
          std::chrono::duration<double>(t1 - t0).count();

        ++nCmod1;
      }

      // --------------------------------------------------------
      // cdiv
      // --------------------------------------------------------

      {
        const auto t0 = clock::now();

        cdiv(
          L,
          j,
          colpointers
        );

        const auto t1 = clock::now();

        tCdiv +=
          std::chrono::duration<double>(t1 - t0).count();

        ++nCdiv;
      }
    }

    // ----------------------------------------------------------
    // Schedule J's own update.
    // ----------------------------------------------------------

    const int width =
      j1 - j0;

    const int len =
      rowpointers[J + 1] -
      rowpointers[J];

    LENGTH[J] =
      len - width;

    if (LENGTH[J] > 0)
    {
      const int next_row =
        rowindices[rowpointers[J] + width];

      const int nextJ =
        SNODE[next_row];

      LINK[J] =
        LINK[nextJ];

      LINK[nextJ] =
        J;
    }
    else
    {
      LENGTH[J] = 0;
      LINK[J] = -1;
    }
  }

  // ------------------------------------------------------------
  // Total time
  // ------------------------------------------------------------

  const auto tTotal1 = clock::now();

  const double tTotal =
    std::chrono::duration<double>(
      tTotal1 - tTotal0
    ).count();

  const double tMeasured =
    tIndMap +
    tCmod2 +
    tCmod1 +
    tCdiv;

  const double tOther =
    tTotal - tMeasured;

  // ------------------------------------------------------------
  // Print timing
  // ------------------------------------------------------------

  Rcpp::Rcout
  << "\nCholesky timing\n"
  << "-------------------------\n"
  << "Total          : " << tTotal << " s\n"
  << "makeIndMap     : " << tIndMap << " s\n"
  << "cmod2          : " << tCmod2 << " s  ("
  << nCmod2 << " calls)\n"
  << "cmod1          : " << tCmod1 << " s  ("
  << nCmod1 << " calls)\n"
  << "cdiv           : " << tCdiv << " s  ("
  << nCdiv << " calls)\n"
  << "other          : " << tOther << " s\n"
  << std::endl;
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

