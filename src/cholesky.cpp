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


void cmod2_sup(
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

  double* tp = t.begin();
  double* l  = L.begin();

  const int eK = rowpointers[K + 1];

  const int sCol = supernodes[K];
  const int eCol = supernodes[K + 1];

  // ------------------------------------------------------------
  // Process all target columns j of J receiving an update
  // from source supernode K.
  //
  // khead, klen and ncolup are supplied by the scheduler, so
  // there is no need for find_targets().
  // ------------------------------------------------------------

  for (int p = 0; p < ncolup; ++p)
  {
    const int j = rowindices[khead + p];
    const int sz = klen - p;

    // Initialise t.
    for (int i = 0; i < sz; ++i)
      tp[i] = 0.0;

    // ----------------------------------------------------------
    // Contribution from all columns of source supernode K.
    //
    // The last sz entries of each column contain the active
    // suffix needed for target column j.
    // ----------------------------------------------------------

    for (int k = sCol; k < eCol; ++k)
    {
      const int jk = colpointers[k + 1] - sz;
      int ik = jk;

      const double Ljk = l[jk];

      for (int i = sz - 1; i >= 0; --i)
      {
        tp[i] += l[ik++] * Ljk;
      }
    }

    // ----------------------------------------------------------
    // Scatter back into target column j.
    // ----------------------------------------------------------

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
  NumericVector t(N);

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
      cmod2_sup(L, J, K, khead, klen, ncolup, t, indmap,
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


