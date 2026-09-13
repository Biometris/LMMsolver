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

//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <set>
#include <vector>
#include "AuxFun.h"
#include "SparseMatrix.h"
#include "cholesky.h"
#include <R_ext/BLAS.h>

#include <chrono>

#include "matvec_test.h"

using namespace std::chrono;

using namespace Rcpp;
using namespace std;


//#include <Rcpp.h>
//#include <chrono>

using namespace Rcpp;
using namespace std::chrono;

// [[Rcpp::export]]
double axpy_test_large(
    int n = 492,
    int nvec = 40000)
{
  NumericVector t(n, 0.0);
  NumericVector src(n * nvec, 1.0);

  double* tp = t.begin();
  const double* sp = src.begin();

  const double alpha = 2.0;

  auto start = std::chrono::high_resolution_clock::now();

  for (int r = 0; r < nvec; ++r)
  {
    const double* x = sp + r * n;

    for (int i = 0; i < n; ++i)
      tp[i] += x[i] * alpha;
  }

  auto end = std::chrono::high_resolution_clock::now();

  return std::chrono::duration<double>(
    end - start).count();
}

// [[Rcpp::export]]
double axpy_test(int n = 492, int reps = 100000)
{
  NumericVector t(n, 0.0);
  NumericVector src(n);

  // Initialize source
  for (int i = 0; i < n; ++i)
    src[i] = 1.0;

  double alpha = 2.0;

  double* tp = t.begin();
  const double* sp = src.begin();

  auto start = high_resolution_clock::now();

  for (int rep = 0; rep < reps; ++rep)
  {
    for (int i = 0; i < n; ++i)
      tp[i] += sp[i] * alpha;
  }

  auto end = high_resolution_clock::now();

  const double elapsed =
    duration<double>(end - start).count();

  Rcout << "n = " << n
        << " reps = " << reps
        << " total = " << elapsed
        << " s"
        << " per AXPY = " << elapsed / reps
        << " s"
        << std::endl;

  // Return something dependent on the calculation so that
  // the result cannot simply be discarded.
  return t[0];
}

// [[Rcpp::export]]
double axpy_test4(int n = 492, int reps = 100000)
{
  NumericVector t(n, 0.0);
  NumericVector src(n);

  for (int i = 0; i < n; ++i)
    src[i] = 1.0;

  const double alpha = 2.0;

  double* tp = t.begin();
  const double* sp = src.begin();

  auto start = std::chrono::high_resolution_clock::now();

  for (int rep = 0; rep < reps; ++rep)
  {
    int i = 0;

    for (; i + 3 < n; i += 4)
    {
      tp[i    ] += sp[i    ] * alpha;
      tp[i + 1] += sp[i + 1] * alpha;
      tp[i + 2] += sp[i + 2] * alpha;
      tp[i + 3] += sp[i + 3] * alpha;
    }

    for (; i < n; ++i)
      tp[i] += sp[i] * alpha;
  }

  auto end = std::chrono::high_resolution_clock::now();

  const double elapsed =
    std::chrono::duration<double>(end - start).count();

  Rcout << "total = " << elapsed
        << " s, per AXPY = " << elapsed / reps
        << " s\n";

  return t[0];
}




// [[Rcpp::export]]
double test_dgemm(int m = 492, int k = 246, int n = 492)
{
  arma::mat A = arma::randu<arma::mat>(m, k);
  arma::mat B = arma::randu<arma::mat>(n, k);
  arma::mat U(m, n, arma::fill::zeros);

  const double alpha = 1.0;
  const double beta  = 0.0;

  const int lda = m;
  const int ldb = n;
  const int ldu = m;

  auto t0 = high_resolution_clock::now();

  for (int rep = 0; rep < 20; ++rep)
  {
    F77_CALL(dgemm)(
        "N", "T",
        &m, &n, &k,
        &alpha,
        A.memptr(), &lda,
        B.memptr(), &ldb,
        &beta,
        U.memptr(), &ldu
    FCONE FCONE
    );
  }

  auto t1 = high_resolution_clock::now();

  return duration<double>(t1 - t0).count() / 20.0;
}


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

// [[Rcpp::export]]
void arma_info()
{
  Rcpp::Rcout << "Armadillo version: "
              << arma::arma_version::as_string() << "\n";

  Rcpp::Rcout << "BLAS enabled: "
              << arma::arma_config::blas << "\n";

  Rcpp::Rcout << "LAPACK enabled: "
              << arma::arma_config::lapack << "\n";

  Rcpp::Rcout << "BLAS wrapper: "
              << arma::arma_config::wrapper << "\n";
}

// j is current column in Supernode J
void cmod1(NumericVector& L, int j, int J,
           const IntegerVector& supernodes,
           const IntegerVector& colpointers)
{
  const int& s = colpointers[j];
  const int& e = colpointers[j+1];
  // for all columns in supernode J left to j:
  for (int k=supernodes[J];k<j;k++)
  {
    const int& jk = colpointers[k] + (j-k);
    int ik = jk;
    const double& Ljk = L[jk];
    for (int ij=s; ij<e; ij++)
    {
       L[ij] -= L[ik++]*Ljk;
       //ik++;
    }
  }
}

// Adjust column j for all columns in supernode K:
inline void cmod2(NumericVector& L, int j, int K, int sz,
           NumericVector& t,
           const IntegerVector& indmap,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  // init t:
  for (int i=0;i<sz;i++)
  {
    t[i] = 0.0;
  }

  const int& sCol = supernodes[K];
  const int& eCol = supernodes[K+1];

  for (int k=sCol; k<eCol; k++)
  {
    int jk = colpointers[k+1]-sz;
    int ik = jk;
    const double& Ljk = L[jk];
    for (int i=sz-1;i>=0;i--)
    {
      t[i] += L[ik++]*Ljk;
      //ik++;
    }
  }

  int r = rowpointers[K+1]-1;
  int ref_pos = colpointers[j+1] - 1;
  for (int i=0;i<sz;i++)
  {
    int ndx = rowindices[r--];
    int pos = ref_pos - indmap[ndx];
    L[pos] -= t[i];
  }
}

void cdiv(NumericVector& L, int j, const IntegerVector& colpointers)
{
  const int& s = colpointers[j];
  const int& e = colpointers[j+1];

  // pivot:
  L[s] = sqrt(L[s]);
  // update column j:
  double Ls = L[s];
  for (int i = s + 1; i < e; i++)
  {
    L[i] /= Ls;
  }
}

// [[Rcpp::export]]
void chol_last_supernode(
    NumericVector& L,
    const IntegerVector& supernodes,
    const IntegerVector& colpointers)
{
  const int J = supernodes.size()-2;
  //Rcout << "J = " << J << endl;

  const int j0 = supernodes[J];
  const int j1 = supernodes[J + 1];

  const int n = j1 - j0;       // should be 492

  // ------------------------------------------------------------
  // 1. Sparse -> dense
  // ------------------------------------------------------------
  auto s0 = high_resolution_clock::now();

  arma::mat A(n, n, arma::fill::zeros);

  double* l = L.begin();

  for (int j = 0; j < n; ++j)
  {
    const int col = j0 + j;
    const int s = colpointers[col];

    for (int i = j; i < n; ++i)
      A(i, j) = A(j, i) = l[s + (i - j)];
  }

  auto s1 = high_resolution_clock::now();

  double time_s2d =
    duration<double>(s1 - s0).count();

  // ------------------------------------------------------------
  // 2. Full dense Cholesky
  // ------------------------------------------------------------
  auto s2 = high_resolution_clock::now();

  A = arma::chol(A).t();

  auto s3 = high_resolution_clock::now();

  double time_chol = duration<double>(s3 - s2).count();


  // ------------------------------------------------------------
  // 3. Dense -> sparse
  // ------------------------------------------------------------
  auto s4 = high_resolution_clock::now();

  for (int j = 0; j < n; ++j)
  {
    const int col = j0 + j;
    const int s = colpointers[col];

    for (int i = j; i < n; ++i)
      l[s + (i - j)] = A(i, j);
  }

  auto s5 = high_resolution_clock::now();

  double time_d2s =
    duration<double>(s5 - s4).count();


  // ------------------------------------------------------------
  // Total
  // ------------------------------------------------------------
  double time_total = duration<double>(s5 - s0).count();

  return;
}


// Update supernode J using supernode K:
//
//     C_J <- C_J - L[I,K] L[J,K]^T
//
// where I is the intersection of the row structures of K and J.
//
// This is an Armadillo prototype. The sparse indexing is done
// outside the matrix multiplication.
//
// K < J
void cmod2_sup(
    NumericVector& L,
    int K,
    int J,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  // ------------------------------------------------------------
  // 1. Find common rows of K and J
  // ------------------------------------------------------------

  auto start1 = high_resolution_clock::now();


  const int sK = rowpointers[K];
  const int eK = rowpointers[K + 1];

  const int sJ = rowpointers[J];
  const int eJ = rowpointers[J + 1];

  std::vector<int> rowsK;
  std::vector<int> rowsJ;

  int k = sK;
  int j = sJ;

  while (k < eK && j < eJ)
  {
    const int rowK = rowindices[k];
    const int rowJ = rowindices[j];

    if (rowK == rowJ)
    {
      rowsK.push_back(k - sK);
      rowsJ.push_back(j - sJ);

      ++k;
      ++j;
    }
    else if (rowK < rowJ)
    {
      ++k;
    }
    else
    {
      ++j;
    }
  }
  auto end1 = high_resolution_clock::now();


  const int h = rowsK.size();
  const int wK = supernodes[K + 1] - supernodes[K];
  const int wJ = supernodes[J + 1] - supernodes[J];

  if (h == 0)
    return;

  // ------------------------------------------------------------
  // 2. Construct
  //
  //     A = L[I,K]       h x wK
  //
  //     B = L[J,K]       wJ x wK
  //
  // For the target supernode, its first wJ row positions are
  // its own columns.
  // ------------------------------------------------------------

  auto start2 = high_resolution_clock::now();

  arma::mat  A(h, wK, arma::fill::zeros);
  arma::mat  B(wJ, wK, arma::fill::zeros);

  double* l = L.begin();

  const int k0 = supernodes[K];
  const int j0 = supernodes[J];

  // ---- A = L[I,K] -------------------------------------------

  for (int c = 0; c < wK; ++c)
  {
    const int col = k0 + c;
    const int s = colpointers[col];

    for (int q = 0; q < h; ++q)
    {
      const int rk = rowsK[q];

      // Column c starts at its diagonal row c.
      if (rk >= c)
        A(q, c) = l[s + rk - c];
    }
  }

  // ---- B = L[J,K] -------------------------------------------
  //
  // Rows of J are its supernode columns:
  //
  //     j0, j0+1, ..., j0+wJ-1
  //
  // Find their positions in K's row structure through the
  // common-row correspondence.
  //
  for (int p = 0; p < wJ; ++p)
  {
    const int rJ = p;  // local row position in J

    // Find this row among the common rows.
    for (int q = 0; q < h; ++q)
    {
      if (rowsJ[q] == rJ)
      {
        const int rk = rowsK[q];

        for (int c = 0; c < wK; ++c)
        {
          const int col = k0 + c;
          const int s = colpointers[col];

          if (rk >= c)
            B(p, c) = l[s + rk - c];
        }

        break;
      }
    }
  }

  auto end2 = high_resolution_clock::now();


  // ------------------------------------------------------------
  // 3. Dense sup-sup update
  //
  //     U = A B^T
  // ------------------------------------------------------------

  auto start3a = high_resolution_clock::now();

  arma::mat Bt = B.t();

  //arma::mat U = A * B.t();

  auto end3a = high_resolution_clock::now();

  auto start3b = high_resolution_clock::now();

  arma::mat U = A * Bt;

  auto end3b = high_resolution_clock::now();


  // ------------------------------------------------------------
  // 4. Scatter U back into J
  //
  // U(q,p) corresponds to
  //
  //   row    = rowindices[sJ + rowsJ[q]]
  //   column = j0 + p
  //
  // Only entries in the lower triangle are stored.
  // ------------------------------------------------------------

  //const int j0 = supernodes[J];
  auto start4 = high_resolution_clock::now();

  for (int q = 0; q < h; ++q)
  {
    const int row =
      rowindices[sJ + rowsJ[q]];

    for (int p = 0; p < wJ; ++p)
    {
      const int col = j0 + p;

      // Only lower triangle
      if (row < col)
        continue;

      const int pos =
        colpointers[col] + (row - col);

      l[pos] -= U(q, p);
    }
  }
  auto end4 = high_resolution_clock::now();
  double elapsed1 = duration<double>(end1 - start1).count();
  double elapsed2 = duration<double>(end2 - start2).count();
  double elapsed3a = duration<double>(end3a - start3a).count();
  double elapsed3b = duration<double>(end3b - start3b).count();
  double elapsed4 = duration<double>(end4 - start4).count();
  Rcout << "K = " << setw(3) << K
        << " J = " << setw(3) << J
        << " struct    " << setw(3) << elapsed1
        << " load      " << setw(3) << elapsed2
        << " tr   " << setw(3) << elapsed3a
        << " matvec_tile " << setw(3) << elapsed3b
        << " write  " << setw(3) << elapsed4
        << std::endl;


  // ------------------------------------------------------------
  // 4. Scatter U back into J
  //
  // U(q,p) corresponds to:
  //
  //     row = rowsJ[q]
  //     col = p
  //
  // Only lower triangle is stored.
  // ------------------------------------------------------------

  //for (int q = 0; q < h; ++q)
  //{
  //  const int rowJ = rowsJ[q];
  //
  //  for (int p = 0; p < wJ; ++p)
  //  {
  //    if (rowJ < p)
  //      continue;

  //    const int col = j0 + p;
  //    const int s = colpointers[col];

      // In column p, local row offset is rowJ - p.
  //    l[s + rowJ - p] -= U(q, p);
  //  }
  //}
  //Rcout << "U = " << setw(12) << U << endl;
}

// K-first prototype:
// Process all columns j of supernode J that receive a contribution from K.
//
// The numerical calculation itself is unchanged: we call the existing
// cmod2() for each target column j.  The only change is the order in
// which the updates are performed.
//
// indmap must already have been constructed for supernode J.

// Process one source supernode K for all columns of target
// supernode J that K contributes to.
//
// This is only a change in the ORDER of the cmod2 calls.
// The numerical operation itself is exactly the original cmod2().

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
  const int sK = rowpointers[K];
  const int eK = rowpointers[K + 1];

  const int k0 = supernodes[K];
  const int k1 = supernodes[K + 1];

  const int j0 = supernodes[J];
  const int j1 = supernodes[J + 1];

  const int wK = k1 - k0;

  // ------------------------------------------------------------
  // Find all target columns j of J to which K contributes.
  // ------------------------------------------------------------
  auto s1 = high_resolution_clock::now();

  std::vector<int> target_j;
  std::vector<int> target_k;

  for (int k = sK; k < eK; ++k)
  {
    const int row = rowindices[k];

    if (row < j0)
      continue;

    if (row >= j1)
      break;

    target_j.push_back(row);
    target_k.push_back(k);
  }
  auto e1 = high_resolution_clock::now();

  const int q = target_j.size();

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
  // of K column k0 + pc.
  // ------------------------------------------------------------
  auto s2 = high_resolution_clock::now();


  std::vector<double> P(wK * maxsz);

  for (int pc = 0; pc < wK; ++pc)
  {
    const int col = k0 + pc;

    const int start =
      colpointers[col + 1] - maxsz;

    for (int r = 0; r < maxsz; ++r)
      P[pc * maxsz + r] = L[start + r];
  }
  auto e2 = high_resolution_clock::now();

  // ------------------------------------------------------------
  // Reuse packed K for all target columns j.
  // ------------------------------------------------------------
  auto s3 = high_resolution_clock::now();

  for (int p = 0; p < q; ++p)
  {
    const int j = target_j[p];

    // Position in K's row structure corresponding to j.
    const int khead = target_k[p];

    // Exactly the same sz as original cmod2().
    const int sz = eK - khead;

    // Offset into the packed suffix.
    const int offset = khead - firstK;

    // ----------------------------------------------------------
    // t = L[I,K] L[j,K]^T
    // ----------------------------------------------------------

    for (int i = 0; i < sz; ++i)
      t[i] = 0.0;

    for (int pc = 0; pc < wK; ++pc)
    {
      const double* Pcol =
        P.data() + pc * maxsz + offset;

      // Pcol[0] corresponds to L[j,k].
      const double Ljk = Pcol[0];


      // Preserve original cmod2() ordering.
      for (int i = sz - 1, r = 0;
           i >= 0;
           --i, ++r)
      {
        t[i] += Pcol[r] * Ljk;
      }
    }

    // ----------------------------------------------------------
    // Scatter exactly as in original cmod2().
    // ----------------------------------------------------------

    int r = eK - 1;

    const int ref_pos =
      colpointers[j + 1] - 1;

    for (int i = 0; i < sz; ++i)
    {
      const int ndx = rowindices[r--];

      const int pos =
        ref_pos - indmap[ndx];

      L[pos] -= t[i];
    }
  }
  auto e3 = high_resolution_clock::now();
  double elapsed1 = duration<double>(e1 - s1).count();
  double elapsed2 = duration<double>(e2 - s2).count();
  double elapsed3 = duration<double>(e3 - s3).count();
  Rcout << "K =  " << setw(3) << K
        << "structure "  << elapsed1
        << "pack "       << elapsed2
        << "calculate "  << elapsed3 << endl;
}

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
  std::vector<int> target_k;

  for (int k = sK; k < eK; ++k)
  {
    const int row = rowindices[k];

    if (row < j0)
      continue;

    if (row >= j1)
      break;

    target_k.push_back(k);
  }

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


void cmod2_Kpanel(
    NumericVector& L,
    int K,
    int J,
    int panelSize,
    NumericVector& t,
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

  // ------------------------------------------------------------
  // Find all target columns j of J to which K contributes.
  //
  // target_j[p] = target column j
  // target_k[p] = position of j in K's row structure
  // ------------------------------------------------------------

  std::vector<int> target_j;
  std::vector<int> target_k;

  for (int k = sK; k < eK; ++k)
  {
    const int row = rowindices[k];

    if (row < j0)
      continue;

    if (row >= j1)
      break;

    target_j.push_back(row);
    target_k.push_back(k);
  }

  const int q = target_j.size();

  if (q == 0)
    return;

  // ------------------------------------------------------------
  // Largest suffix needed.
  // ------------------------------------------------------------

  const int firstK = target_k[0];
  const int maxsz = eK - firstK;

  if (t.size() < maxsz)
    stop("t is too small");

  // ------------------------------------------------------------
  // Process K panel by panel.
  // ------------------------------------------------------------

  for (int panelStart = 0;
       panelStart < wK;
       panelStart += panelSize)
  {
    const int panelEnd =
      std::min(wK, panelStart + panelSize);

    const int panelWidth = panelEnd - panelStart;

    // ----------------------------------------------------------
    // Pack the last maxsz entries of each K column.
    // ----------------------------------------------------------

    std::vector<double> P(panelWidth * maxsz);

    for (int pc = 0; pc < panelWidth; ++pc)
    {
      const int col = k0 + panelStart + pc;

      const int start =
        colpointers[col + 1] - maxsz;

      for (int r = 0; r < maxsz; ++r)
        P[pc * maxsz + r] = L[start + r];
    }

    // ----------------------------------------------------------
    // Use this K-panel for every target column j.
    // ----------------------------------------------------------

    for (int p = 0; p < q; ++p)
    {
      const int j = target_j[p];

      // Equivalent to the original colhead[K].
      const int khead = target_k[p];

      // Exactly the same sz as original cmod2().
      const int sz = eK - khead;

      // The target column's jk lies this far into the packed
      // suffix.
      const int offset = khead - firstK;

      // --------------------------------------------------------
      // t = contribution from this K-panel
      // --------------------------------------------------------

      for (int i = 0; i < sz; ++i)
        t[i] = 0.0;

      for (int pc = 0; pc < panelWidth; ++pc)
      {
        const double* Pcol =
          P.data() + pc * maxsz + offset;

        // Pcol[0] is exactly the original L[jk],
        // hence L[j,k].
        const double Ljk = Pcol[0];

        // IMPORTANT:
        //
        // Original:
        //
        //   i = sz-1 ... 0
        //   t[i] += L[jk + (sz-1-i)] * Ljk
        //
        // Therefore the packed Pcol is traversed forward,
        // while t is traversed backwards.
        //
        for (int i = sz - 1, r = 0;
             i >= 0;
             --i, ++r)
        {
          t[i] += Pcol[r] * Ljk;
        }
      }

      // --------------------------------------------------------
      // Scatter exactly as in original cmod2().
      // --------------------------------------------------------

      int r = eK - 1;

      const int ref_pos =
        colpointers[j + 1] - 1;

      for (int i = 0; i < sz; ++i)
      {
        const int ndx = rowindices[r--];

        const int pos =
          ref_pos - indmap[ndx];

        L[pos] -= t[i];
      }
    }
  }
}

void cmod2_Kpanel_blocked(
    NumericVector& L,
    int K,
    int J,
    int panelSize,
    int jBlockSize,
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

  // ------------------------------------------------------------
  // Find the target columns j of J affected by K.
  //
  // target_j[p] = actual column j
  // target_k[p] = corresponding position in K's row structure
  //
  // This target_k is exactly the position that would be
  // colhead[K] in the original j-first algorithm.
  // ------------------------------------------------------------

  std::vector<int> target_j;
  std::vector<int> target_k;

  for (int k = sK; k < eK; ++k)
  {
    const int r = rowindices[k];

    if (r < j0)
      continue;

    if (r >= j1)
      break;

    target_j.push_back(r);
    target_k.push_back(k);
  }

  const int q = target_j.size();

  if (q == 0)
    return;

  // ------------------------------------------------------------
  // Process K in panels.
  // ------------------------------------------------------------

  for (int kp = 0; kp < wK; kp += panelSize)
  {
    const int kpEnd =
      std::min(wK, kp + panelSize);

    // ----------------------------------------------------------
    // Process several target columns of J at once.
    // ----------------------------------------------------------

    for (int jp = 0; jp < q; jp += jBlockSize)
    {
      const int jpEnd =
        std::min(q, jp + jBlockSize);

      const int jb = jpEnd - jp;

      // --------------------------------------------------------
      // Maximum active size in this J block.
      //
      // target_k increases with j, so earlier target columns
      // generally have the largest sz.
      // --------------------------------------------------------

      int maxsz = 0;

      for (int p = jp; p < jpEnd; ++p)
      {
        const int sz = eK - target_k[p];

        if (sz > maxsz)
          maxsz = sz;
      }

      // --------------------------------------------------------
      // T is a small dense block:
      //
      //   jb target columns x maxsz active rows
      //
      // T[p * maxsz + i]
      //
      // contains the contribution of this K-panel to target
      // column target_j[p].
      // --------------------------------------------------------

      std::vector<double> T(jb * maxsz, 0.0);

      // --------------------------------------------------------
      // Block multiplication:
      //
      //   K-panel -> J-block -> K-column -> row
      //
      // --------------------------------------------------------

      for (int kc = kp; kc < kpEnd; ++kc)
      {
        const int col = k0 + kc;

        const int colEnd = colpointers[col + 1];

        for (int p = jp; p < jpEnd; ++p)
        {
          const int localP = p - jp;

          const int sz =
            eK - target_k[p];

          // This is exactly the original
          //
          //   jk = colpointers[col+1] - sz
          //
          // from cmod2().
          const int jk =
            colEnd - sz;

          const double Ljk = L[jk];

          double* tp =
            T.data() + localP * maxsz;

          // ----------------------------------------------------
          // Preserve exactly the original cmod2() ordering:
          //
          //   i = sz-1 ... 0
          //
          //   t[i] += L[jk + r] * Ljk
          //
          // where r runs forward.
          // ----------------------------------------------------

          int r = 0;

          for (int i = sz - 1;
               i >= 0;
               --i, ++r)
          {
            tp[i] += L[jk + r] * Ljk;
          }
        }
      }

      // --------------------------------------------------------
      // Scatter this J-block immediately.
      // --------------------------------------------------------

      for (int p = jp; p < jpEnd; ++p)
      {
        const int localP = p - jp;

        const int j =
          target_j[p];

        const int sz =
          eK - target_k[p];

        const double* tp =
          T.data() + localP * maxsz;

        int r = eK - 1;

        const int ref_pos =
          colpointers[j + 1] - 1;

        for (int i = 0; i < sz; ++i)
        {
          const int ndx =
            rowindices[r--];

          const int pos =
            ref_pos - indmap[ndx];

          L[pos] -= tp[i];
        }
      }
    }
  }
}



void cmod2_4(
    NumericVector& L,
    int K,
    int J,
    NumericVector& t1,
    NumericVector& t2,
    NumericVector& t3,
    NumericVector& t4,
    const IntegerVector& indmap,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  const int sK = rowpointers[K];
  const int eK = rowpointers[K + 1];

  const int j0 = supernodes[J];
  const int j1 = supernodes[J + 1];

  const int k0 = supernodes[K];
  const int k1 = supernodes[K + 1];

  // ------------------------------------------------------------
  // Find target columns of J affected by K.
  // ------------------------------------------------------------

  std::vector<int> target_j;
  std::vector<int> target_k;

  for (int k = sK; k < eK; ++k)
  {
    const int row = rowindices[k];

    if (row < j0)
      continue;

    if (row >= j1)
      break;

    target_j.push_back(row);
    target_k.push_back(k);
  }

  const int q = target_j.size();

  // ------------------------------------------------------------
  // Process four consecutive target columns at a time.
  // ------------------------------------------------------------

  for (int p = 0; p < q; p += 4)
  {
    const int n = std::min(4, q - p);

    // For this optimized kernel we require consecutive columns.
    // If not, fall back to the original cmod2().
    bool consecutive = true;

    for (int u = 1; u < n; ++u)
    {
      if (target_j[p + u] != target_j[p] + u)
      {
        consecutive = false;
        break;
      }
    }

    if (!consecutive || n < 4)
    {
      for (int u = 0; u < n; ++u)
      {
        const int j = target_j[p + u];
        const int sz = eK - target_k[p + u];

        for (int i = 0; i < sz; ++i)
          t1[i] = 0.0;

        for (int col = k0; col < k1; ++col)
        {
          const int jk = colpointers[col + 1] - sz;
          const double Ljk = L[jk];

          int ik = jk;

          for (int i = sz - 1; i >= 0; --i)
            t1[i] += L[ik++] * Ljk;
        }

        int r = eK - 1;
        const int ref_pos = colpointers[j + 1] - 1;

        for (int i = 0; i < sz; ++i)
        {
          const int ndx = rowindices[r--];
          const int pos = ref_pos - indmap[ndx];

          L[pos] -= t1[i];
        }
      }

      continue;
    }

    // ----------------------------------------------------------
    // Four consecutive columns.
    // ----------------------------------------------------------

    const int jA = target_j[p];

    const int headA = target_k[p];
    const int szA = eK - headA;

    const int szB = szA - 1;
    const int szC = szA - 2;
    const int szD = szA - 3;

    // ----------------------------------------------------------
    // Initialize the four accumulators.
    // ----------------------------------------------------------

    for (int i = 0; i < szA; ++i)
      t1[i] = 0.0;

    for (int i = 0; i < szB; ++i)
      t2[i] = 0.0;

    for (int i = 0; i < szC; ++i)
      t3[i] = 0.0;

    for (int i = 0; i < szD; ++i)
      t4[i] = 0.0;

    // ----------------------------------------------------------
    // Main calculation.
    //
    // For one K-column:
    //
    //   jA starts at x[0]
    //   jA+1 starts at x[1]
    //   jA+2 starts at x[2]
    //   jA+3 starts at x[3]
    //
    // The four source vectors therefore overlap.
    //
    // We exploit that overlap explicitly.
    // ----------------------------------------------------------

    for (int col = k0; col < k1; ++col)
    {
      // Start of the suffix for jA.
      const int jk = colpointers[col + 1] - szA;

      // Values corresponding to the four target rows.
      const double LjA = L[jk];
      const double LjB = L[jk + 1];
      const double LjC = L[jk + 2];
      const double LjD = L[jk + 3];

      // --------------------------------------------------------
      // Common part.
      //
      // x[r] is reused for the four updates.
      // --------------------------------------------------------

      const int common = szD;

      for (int r = 0; r < common; ++r)
      {
        const double x0 = L[jk + r];
        const double x1 = L[jk + r + 1];
        const double x2 = L[jk + r + 2];
        const double x3 = L[jk + r + 3];

        t1[szA - 1 - r] += x0 * LjA;
        t2[szB - 1 - r] += x1 * LjB;
        t3[szC - 1 - r] += x2 * LjC;
        t4[szD - 1 - r] += x3 * LjD;
      }

      // --------------------------------------------------------
      // Remaining tails.
      // --------------------------------------------------------

      for (int r = common; r < szC; ++r)
      {
        const double x0 = L[jk + r];
        const double x1 = L[jk + r + 1];
        const double x2 = L[jk + r + 2];

        t1[szA - 1 - r] += x0 * LjA;
        t2[szB - 1 - r] += x1 * LjB;
        t3[szC - 1 - r] += x2 * LjC;
      }

      for (int r = szC; r < szB; ++r)
      {
        const double x0 = L[jk + r];
        const double x1 = L[jk + r + 1];

        t1[szA - 1 - r] += x0 * LjA;
        t2[szB - 1 - r] += x1 * LjB;
      }

      for (int r = szB; r < szA; ++r)
      {
        t1[szA - 1 - r] += L[jk + r] * LjA;
      }
    }

    // ----------------------------------------------------------
    // Write the four columns back.
    // ----------------------------------------------------------

    // ---- jA --------------------------------------------------

    {
      int r = eK - 1;
      const int ref_pos = colpointers[jA + 1] - 1;

      for (int i = 0; i < szA; ++i)
      {
        const int ndx = rowindices[r--];
        const int pos = ref_pos - indmap[ndx];

        L[pos] -= t1[i];
      }
    }

    // ---- jA + 1 ----------------------------------------------

    {
      const int j = jA + 1;

      int r = eK - 1;
      const int ref_pos = colpointers[j + 1] - 1;

      for (int i = 0; i < szB; ++i)
      {
        const int ndx = rowindices[r--];
        const int pos = ref_pos - indmap[ndx];

        L[pos] -= t2[i];
      }
    }

    // ---- jA + 2 ----------------------------------------------

    {
      const int j = jA + 2;

      int r = eK - 1;
      const int ref_pos = colpointers[j + 1] - 1;

      for (int i = 0; i < szC; ++i)
      {
        const int ndx = rowindices[r--];
        const int pos = ref_pos - indmap[ndx];

        L[pos] -= t3[i];
      }
    }

    // ---- jA + 3 ----------------------------------------------

    {
      const int j = jA + 3;

      int r = eK - 1;
      const int ref_pos = colpointers[j + 1] - 1;

      for (int i = 0; i < szD; ++i)
      {
        const int ndx = rowindices[r--];
        const int pos = ref_pos - indmap[ndx];

        L[pos] -= t4[i];
      }
    }
  }
}


void cholesky(NumericVector& L,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  auto start_chol = high_resolution_clock::now();

  const int N = colpointers.size() - 1;
  const int Nsupernodes = supernodes.size()-1;

  // linked lists, see section 4.2 Ng and Peyton
  IntegerVector HEAD(N,-1);
  IntegerVector LINK(Nsupernodes,-1);

  IntegerVector colhead = clone(rowpointers);
  for (int J=0; J<Nsupernodes;J++)
  {
    int szNode = supernodes[J+1] - supernodes[J];
    colhead[J] += szNode-1;
    if (colhead[J] < rowpointers[J+1]-1)
    {
      int rNdx = rowindices[colhead[J]+1];
      insert(HEAD, LINK, rNdx, J);
    }
  }


  IntegerVector indmap(N,0);
  NumericVector t(N);
  int last_sn = Nsupernodes-1;

  // for each supernode J
  for (int J=0; J<Nsupernodes;J++) {
    // used to test new implementation:
    // last supernode..
    NumericVector copy_L = clone(L);

    if (J == Nsupernodes-1) {

      Rcout << "Left looking for last supernode:..." << endl;
      auto start_cmod2_sup = high_resolution_clock::now();

      // Phase 1
      makeIndMap(indmap, J, rowpointers, rowindices);

      for (int K = 0; K < J; ++K)
      {
        cmod2_Kfirst4(copy_L, K, J, indmap, supernodes, rowpointers, colpointers, rowindices);
      }

      auto end_cmod2_sup = high_resolution_clock::now();
      double elapsed_cmod2_sup = duration<double>(end_cmod2_sup - start_cmod2_sup).count();

      //  cmod2_sup(copy_L, K, J, supernodes, rowpointers, colpointers, rowindices);
      //}
      Rcout << "Supernode " << setw(3) << J
            << " cmod2_Kfirst: " << setw(4) << elapsed_cmod2_sup
            << std::endl;
    }

    long long calls = 0;
    long long work = 0;

    int sz_node = supernodes[J+1] - supernodes[J];

    // Phase 1
    makeIndMap(indmap, J, rowpointers, rowindices);
    auto start_cmod2 = high_resolution_clock::now();

    for (int j=supernodes[J];j<supernodes[J+1];j++)
    {
      int K = HEAD[j];

      while (K!=-1)
      {

        int nextK = LINK[K];
        int sz = rowpointers[K+1] - colhead[K];

        const int wK = supernodes[K+1] - supernodes[K];
        calls++;
        work += static_cast<long long>(sz) * wK;
        cmod2(L, j, K, sz, t, indmap, supernodes, rowpointers, colpointers, rowindices);

        colhead[K]++;
        if (colhead[K] < rowpointers[K+1])
        {
          int rNdx = rowindices[colhead[K]];
          insert(HEAD, LINK, rNdx, K);
        }
        K = nextK;
      }
      HEAD[j] = -1;
    }

    auto end_cmod2 = high_resolution_clock::now();
    double elapsed_cmod2 = duration<double>(end_cmod2 - start_cmod2).count();

    // update
    colhead[J]++;

    if (J == last_sn) {
      Rcout << "Comparison using matvec tile vs original code" << endl;

      auto start_chol_last = high_resolution_clock::now();

      chol_last_supernode(copy_L, supernodes, colpointers);

      auto end_chol_last = high_resolution_clock::now();

      double elapsed_chol_last = duration<double>(end_chol_last - start_chol_last).count();

      Rcout << "Supernode " << setw(3) << J
            << " sznode: " << setw(3) << sz_node
            << " cmod2:  " << setw(3) << elapsed_cmod2
            << " cmod_arma: " << setw(4) << elapsed_chol_last
            << std::endl;
    }

    // Phase 2
    auto start = high_resolution_clock::now();
    for (int j=supernodes[J];j<supernodes[J+1];j++) {
      cmod1(L, j, J, supernodes, colpointers);
      cdiv(L, j, colpointers);
    }
    auto end = high_resolution_clock::now();
    double elapsed = duration<double>(end - start).count();
    if (J == last_sn) {
      if (!Rcpp::is_true(Rcpp::all(Rcpp::abs(copy_L - L) < 1e-09)))
        Rcpp::stop("copy_L and L differ");

    }

    if (sz_node > 2) {
      Rcout << "Supernode " << setw(3) << J << " sznode: " << setw(3) << sz_node <<
        " cmod2: " << setw(4) << elapsed_cmod2 << " cmod1/cdiv: " << setw(4) << elapsed << endl;
    }

  }
  auto end_chol = high_resolution_clock::now();
  double elapsed_total = duration<double>(end_chol - start_chol).count();
  Rcout << endl << "Total time Cholesky: " << elapsed_total << endl << endl;
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
