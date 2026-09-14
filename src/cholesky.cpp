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


/*
void cmod2_Kfirst4(
    NumericVector& L,
    int K,
    int J,
    const IntegerVector& indmap,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  const int sK = rowpointers[K];
  const int eK = rowpointers[K + 1];

  const int k0 = supernodes[K];
  const int k1 = supernodes[K + 1];

  const int j0 = supernodes[J];
  const int j1 = supernodes[J + 1];

  const int wK = k1 - k0;
  const int q  = j1 - j0;

  // ------------------------------------------------------------
  // Find target positions in K's row structure.
  // ------------------------------------------------------------
  std::vector<int> target_k = find_targets(K, J, supernodes, rowpointers, rowindices);

  const int nTarget = target_k.size();

  if (nTarget == 0)
    return;

  // ------------------------------------------------------------
  // Check the supernode structure required by the 4-way kernel.
  //
  // Consecutive target columns must correspond to consecutive
  // positions in K's row structure.
  // ------------------------------------------------------------
  for (int p = 0; p + 1 < nTarget; ++p)
  {
    if (target_k[p + 1] != target_k[p] + 1)
      stop("Unexpected target_k spacing in cmod2_Kfirst4");
  }

  // ------------------------------------------------------------
  // Largest K suffix needed.
  // ------------------------------------------------------------
  const int firstK = target_k[0];
  const int maxsz  = eK - firstK;

  const double* lrp = L.begin();
  double* lp = L.begin();

  // ------------------------------------------------------------
  // Pack K once.
  //
  // P[pc * maxsz + r] contains the suffix of column k0 + pc.
  // ------------------------------------------------------------
  std::vector<double> P(wK * maxsz);

  for (int pc = 0; pc < wK; ++pc)
  {
    const int col = k0 + pc;

    const int start =
      colpointers[col + 1] - maxsz;

    double* Pcol =
      P.data() + pc * maxsz;

    for (int r = 0; r < maxsz; ++r)
      Pcol[r] = lrp[start + r];
  }

  // ------------------------------------------------------------
  // Four target accumulators.
  // ------------------------------------------------------------
  std::vector<double> t0(maxsz);
  std::vector<double> t1(maxsz);
  std::vector<double> t2(maxsz);
  std::vector<double> t3(maxsz);

  int p0 = 0;

  // ------------------------------------------------------------
  // Four target columns at a time.
  // ------------------------------------------------------------
  for (; p0 + 3 < nTarget; p0 += 4)
  {
    const int khead0 = target_k[p0];

    const int offset =
      khead0 - firstK;

    const int sz0 = eK - target_k[p0];
    const int sz1 = eK - target_k[p0 + 1];
    const int sz2 = eK - target_k[p0 + 2];
    const int sz3 = eK - target_k[p0 + 3];

    // This should follow from the consecutive target_k positions.
    if (sz1 != sz0 - 1 ||
        sz2 != sz0 - 2 ||
        sz3 != sz0 - 3)
    {
      stop("Unexpected target suffix lengths in cmod2_Kfirst4");
    }

    // ----------------------------------------------------------
    // Zero accumulators.
    // ----------------------------------------------------------
    for (int i = 0; i < sz0; ++i)
      t0[i] = 0.0;

    for (int i = 0; i < sz1; ++i)
      t1[i] = 0.0;

    for (int i = 0; i < sz2; ++i)
      t2[i] = 0.0;

    for (int i = 0; i < sz3; ++i)
      t3[i] = 0.0;

    // ----------------------------------------------------------
    // K-first calculation.
    //
    // For four consecutive target columns, the mapping simplifies
    // so that all four targets use the same P position for a
    // given i:
    //
    //     x = Pcol[maxsz - 1 - i]
    //
    // while their Ljk values are:
    //
    //     Pcol[offset + 0]
    //     Pcol[offset + 1]
    //     Pcol[offset + 2]
    //     Pcol[offset + 3]
    // ----------------------------------------------------------
    for (int pc = 0; pc < wK; ++pc)
    {
      const double* Pcol =
        P.data() + pc * maxsz;

      const double Lj0k = Pcol[offset    ];
      const double Lj1k = Pcol[offset + 1];
      const double Lj2k = Pcol[offset + 2];
      const double Lj3k = Pcol[offset + 3];

      // --------------------------------------------------------
      // Common part: all four targets.
      // --------------------------------------------------------
      for (int i = 0; i < sz3; ++i)
      {
        const double x =
          Pcol[maxsz - 1 - i];

        t0[i] += x * Lj0k;
        t1[i] += x * Lj1k;
        t2[i] += x * Lj2k;
        t3[i] += x * Lj3k;
      }

      // --------------------------------------------------------
      // Remaining part: targets 0,1,2.
      // --------------------------------------------------------
      for (int i = sz3; i < sz2; ++i)
      {
        const double x =
          Pcol[maxsz - 1 - i];

        t0[i] += x * Lj0k;
        t1[i] += x * Lj1k;
        t2[i] += x * Lj2k;
      }

      // --------------------------------------------------------
      // Remaining part: targets 0,1.
      // --------------------------------------------------------
      for (int i = sz2; i < sz1; ++i)
      {
        const double x =
          Pcol[maxsz - 1 - i];

        t0[i] += x * Lj0k;
        t1[i] += x * Lj1k;
      }

      // --------------------------------------------------------
      // Remaining part: target 0.
      // --------------------------------------------------------
      for (int i = sz1; i < sz0; ++i)
      {
        const double x =
          Pcol[maxsz - 1 - i];

        t0[i] += x * Lj0k;
      }
    }

    // ----------------------------------------------------------
    // Scatter target 0.
    // ----------------------------------------------------------
    {
      const int j = j0 + p0;

      int r = eK - 1;
      const int ref_pos = colpointers[j + 1] - 1;

      for (int i = 0; i < sz0; ++i)
      {
        const int ndx = rowindices[r--];
        const int pos = ref_pos - indmap[ndx];

        lp[pos] -= t0[i];
      }
    }

    // ----------------------------------------------------------
    // Scatter target 1.
    // ----------------------------------------------------------
    {
      const int j = j0 + p0 + 1;

      int r = eK - 1;
      const int ref_pos = colpointers[j + 1] - 1;

      for (int i = 0; i < sz1; ++i)
      {
        const int ndx = rowindices[r--];
        const int pos = ref_pos - indmap[ndx];

        lp[pos] -= t1[i];
      }
    }

    // ----------------------------------------------------------
    // Scatter target 2.
    // ----------------------------------------------------------
    {
      const int j = j0 + p0 + 2;

      int r = eK - 1;
      const int ref_pos = colpointers[j + 1] - 1;

      for (int i = 0; i < sz2; ++i)
      {
        const int ndx = rowindices[r--];
        const int pos = ref_pos - indmap[ndx];

        lp[pos] -= t2[i];
      }
    }

    // ----------------------------------------------------------
    // Scatter target 3.
    // ----------------------------------------------------------
    {
      const int j = j0 + p0 + 3;

      int r = eK - 1;
      const int ref_pos = colpointers[j + 1] - 1;

      for (int i = 0; i < sz3; ++i)
      {
        const int ndx = rowindices[r--];
        const int pos = ref_pos - indmap[ndx];

        lp[pos] -= t3[i];
      }
    }
  }

  // ------------------------------------------------------------
  // Remaining 1-3 target columns.
  // Use the original packed K-first calculation.
  // ------------------------------------------------------------
  std::vector<double> t(maxsz);

  for (; p0 < nTarget; ++p0)
  {
    const int j = j0 + p0;

    const int khead = target_k[p0];

    const int sz =
      eK - khead;

    const int offset =
      khead - firstK;

    for (int i = 0; i < sz; ++i)
      t[i] = 0.0;

    for (int pc = 0; pc < wK; ++pc)
    {
      const double* Pcol =
        P.data() + pc * maxsz + offset;

      const double Ljk = Pcol[0];

      for (int i = sz - 1, r = 0;
           i >= 0;
           --i, ++r)
      {
        t[i] += Pcol[r] * Ljk;
      }
    }

    int r = eK - 1;
    const int ref_pos = colpointers[j + 1] - 1;

    for (int i = 0; i < sz; ++i)
    {
      const int ndx = rowindices[r--];
      const int pos = ref_pos - indmap[ndx];

      lp[pos] -= t[i];
    }
  }
}





// ------------------------------------------------------------
// Process one source supernode K for all columns of target
// supernode J that K contributes to.
//
// Same numerical calculation as original cmod2(), but with
// the target columns processed together after packing K.
// ------------------------------------------------------------
void cmod2_Kfirst(
    NumericVector& L,
    int K,
    int J,
    NumericVector& t,
    const IntegerVector& indmap,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  const int eK = rowpointers[K + 1];

  const int k0 = supernodes[K];
  const int k1 = supernodes[K + 1];
  const int wK = k1 - k0; // number of columns in K

  const int j0 = supernodes[J];

  // ------------------------------------------------------------
  // Find all target positions in K's row structure.
  // ------------------------------------------------------------
  std::vector<int> target_k = find_targets(K, J, supernodes, rowpointers, rowindices);

  const int q = target_k.size();

  if (q == 0)
    return;

  // ------------------------------------------------------------
  // Largest suffix needed: first target column.
  // ------------------------------------------------------------
  const int firstK = target_k[0];
  const int maxsz = eK - firstK;

  if (t.size() < maxsz)
    stop("t is too small");

  // ------------------------------------------------------------
  // Pack ALL of K once.
  //
  // P[pc * maxsz + r] contains the last maxsz entries
  // of column k0 + pc.
  // ------------------------------------------------------------
  std::vector<double> P(wK * maxsz);

  const double* lp = L.begin();

  for (int k = 0; k < wK; ++k)
  {
    const int col = k0 + k;

    const int start = colpointers[col + 1] - maxsz;

    double* Pcol = P.data() + k * maxsz;

    for (int r = 0; r < maxsz; ++r)
      Pcol[r] = lp[start + r];
  }

  // ------------------------------------------------------------
  // Reuse packed K for all target columns.
  // ------------------------------------------------------------
  double* tp = t.begin();

  for (int p = 0; p < q; ++p)
  {
    const int j = j0 + p;

    // Position in K's row structure corresponding to j.
    const int khead = target_k[p];

    // Exactly the same sz as original cmod2().
    const int sz = eK - khead;

    // Offset into packed suffix.
    const int offset = khead - firstK;

    // ----------------------------------------------------------
    // t = L[I,K] L[j,K]^T
    // ----------------------------------------------------------

    for (int i = 0; i < sz; ++i) {
      tp[i] = 0.0;
    }

    for (int k = 0; k < wK; ++k)
    {
      const double* Pcol = P.data() + k * maxsz + offset;

      // Pcol[0] corresponds to L[j,k].
      const double Ljk = Pcol[0];

      // Preserve original cmod2() ordering.
      for (int i = sz - 1, r = 0; i >= 0; --i, ++r)
      {
        tp[i] += Pcol[r] * Ljk;
      }
    }

    // ----------------------------------------------------------
    // Scatter exactly as in original cmod2().
    // ----------------------------------------------------------

    int r = eK - 1;

    const int ref_pos = colpointers[j + 1] - 1;

    for (int i = 0; i < sz; ++i)
    {
      const int ndx = rowindices[r--];
      const int pos = ref_pos - indmap[ndx];
      L[pos] -= tp[i];
    }
  }
}
*/



/*
 *
// [[Rcpp::export]]
double logdet(Rcpp::S4 obj, NumericVector lambda)
{
  //Rcpp::S4 obj(arg);
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericVector L = obj.slot("entries");
  NumericMatrix P = obj.slot("P");
  //NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

  // define matrix L (lower triangle matrix values)
  const int sz = P.nrow();
  const int n_prec_mat = P.ncol();
  //NumericVector L(sz, 0.0);
  for (int i=0;i<sz;i++)
  {
    L[i] = 0.0;
  }
  for (int k=0;k<n_prec_mat;k++)
  {
    NumericMatrix::Column Pk = P(_, k);
    double alpha = lambda[k];
    for (int i=0;i<sz;i++)
    {
      L[i] += alpha*Pk[i];
    }
  }
  cholesky(L, supernodes, rowpointers, colpointers, rowindices);
  return logdet(L, colpointers);
}

*/

/*
// Just to show the structure of the sparse cholesky matrix with supernodes.
// [[Rcpp::export]]
NumericMatrix PrintCholesky(Rcpp::S4 obj)
{
  Rcout << "Class: " << as<std::string>(obj.attr("class")) << std::endl;

  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  IntegerVector pivot = obj.slot("pivot");
  IntegerVector invpivot = obj.slot("invpivot");

  NumericVector L = obj.slot("entries");

  const int Nsupernodes = supernodes.size()-1;
  const int N = colpointers.size() - 1;
  NumericMatrix A(N, N);
  for (int J=0; J<Nsupernodes;J++)
  {
    int s = rowpointers[J];
    Rcout << "Supernode: " << J << endl;
    for (int j=supernodes[J]; j<supernodes[J+1]; j++)
    {
      Rcout << "  Column: " << j << endl;
      int k = s;
      for (int ndx = colpointers[j]; ndx < colpointers[j+1]; ndx++)
      {
        int i = rowindices[k++];
        Rcout << "    row: " << i << " (ndx or key " << ndx << ")" << endl;
        A(i, j) = L[ndx];
      }
      s++;
    }
  }
  return A;
}
*/
