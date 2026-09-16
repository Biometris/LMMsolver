#ifndef DIFFCHOLESKYAUXFUN_HEADER
#define DIFFCHOLESKYAUXFUN_HEADER

#include <Rcpp.h>
#include <set>
#include <vector>

using namespace Rcpp;
using namespace std;

// Transform to C++ Notation indices
void transf2C(IntegerVector& ndx);

IntegerVector GetIntVector(Rcpp::S4 obj, const String& slotName, int ArrayIndexing);

NumericVector GetNumericVector(Rcpp::S4 obj, const String& slotName);

void insert(IntegerVector& HEAD, IntegerVector& LINK, int i, int J);

template <int M, int N>
inline void matmul_block(
    const double* A,
    const double* Bt,
    double* C,
    int lda,
    int ldb,
    int ldc,
    int K)
{
  double c[M][N] = {};

  for (int k = 0; k < K; ++k)
  {
    const double* b = Bt + k * ldb;

    #pragma GCC unroll M
    for (int i = 0; i < M; ++i)
    {
      const double a = A[i * lda + k];

      #pragma GCC unroll N
      for (int j = 0; j < N; ++j)
        c[i][j] += a * b[j];
    }
  }

  for (int i = 0; i < M; ++i)
  {
    for (int j = 0; j < N; ++j)
      C[i * ldc + j] = c[i][j];
  }
}


// ============================================================
// Update an M x N block
//
// A  : M x srcWidth
// Bt : srcWidth x N
// C  : M x N
//
// Performs:
//     C = A Bt
//
// and scatters C into the sparse L matrix.
// ============================================================

template <int M, int N>
inline void update_block(
    double* l,
    std::vector<double>& A,
    std::vector<double>& Bt,
    int srcWidth,
    int r0,
    int p0,
    int klen,
    int eK,
    const std::vector<int>& srcBase,
    int khead,
    const IntegerVector& indmap,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  // ----------------------------------------------------------
  // Pack A = X[r0:r0+M-1, :]
  // ----------------------------------------------------------

  for (int k = 0; k < srcWidth; ++k)
  {
    const int base = srcBase[k];

    for (int i = 0; i < M; ++i)
      A[i * srcWidth + k] =
        l[base + r0 + i];
  }

  // ----------------------------------------------------------
  // C = A Bt
  // ----------------------------------------------------------

  double C[M * N];

  matmul_block<M, N>(
      A.data(),
      Bt.data(),
      C,
      srcWidth,
      N,
      N,
      srcWidth
  );

  // ----------------------------------------------------------
  // Scatter valid lower-triangular elements
  // ----------------------------------------------------------

  for (int j = 0; j < N; ++j)
  {
    const int p = p0 + j;

    // On the diagonal block only rows r >= p are valid.
    const int imin =
      (r0 == p0) ? j : 0;

    const int target =
      rowindices[khead + p];

    const int ref_pos =
      colpointers[target + 1] - 1;

    for (int i = imin; i < M; ++i)
    {
      const int r = r0 + i;

      const int q =
        klen - r - 1;

      const int ndx =
        rowindices[eK - 1 - q];

      const int pos =
        ref_pos - indmap[ndx];

      l[pos] -= C[i * N + j];
    }
  }
}

#endif
