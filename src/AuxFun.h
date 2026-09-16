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

// ============================================================
// Reverse/AD matrix multiplication
//
// A  : M x K
// Bt : K x N
// FC : M x N
//
// FA  = FC Bt^T    -> M x K
// FBt = A^T FC     -> K x N
//
// M and N are compile-time constants.
// ============================================================

template <int M, int N>
inline void ADmatmul_block(
    const double* A,
    const double* Bt,
    const double* FC,
    double* FA,
    double* FBt,
    int lda,
    int ldbt,
    int ldfc,
    int ldfa,
    int ldfbt,
    int K)
{
  // ----------------------------------------------------------
  // FA = FC Bt^T
  // ----------------------------------------------------------

  for (int k = 0; k < K; ++k)
  {
    const double* b = Bt + k * ldbt;

    #pragma GCC unroll M
    for (int i = 0; i < M; ++i)
    {
      double sum = 0.0;

      #pragma GCC unroll N
      for (int j = 0; j < N; ++j)
        sum += FC[i * ldfc + j] * b[j];

      FA[i * ldfa + k] = sum;
    }
  }

  // ----------------------------------------------------------
  // FBt = A^T FC
  // ----------------------------------------------------------

  for (int k = 0; k < K; ++k)
  {
    #pragma GCC unroll N
    for (int j = 0; j < N; ++j)
    {
      double sum = 0.0;

      #pragma GCC unroll M
      for (int i = 0; i < M; ++i)
        sum += A[i * lda + k] * FC[i * ldfc + j];

      FBt[k * ldfbt + j] = sum;
    }
  }
}

// ============================================================
// Reverse/AD M x N block update
//
// Forward block:
//     C = A Bt
//
// Reverse:
//     FA  -= FC Bt^T
//     FBt -= A^T FC
//
// q = 0 denotes the bottom entry of the off-diagonal part
// of source supernode K.
//
// A   : M x srcWidth
// Bt  : srcWidth x N
// FC  : M x N
// FA  : M x srcWidth
// FBt : srcWidth x N
// ============================================================

template <int M, int N>
inline void ADupdate_block(
    const double* l,
    double* f,
    std::vector<double>& A,
    std::vector<double>& Bt,
    std::vector<double>& FC,
    std::vector<double>& FA,
    std::vector<double>& FBt,
    int srcWidth,
    int r0,
    int p0,
    int done,
    int eK,
    const std::vector<int>& srcEnd,
    const IntegerVector& indmap,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices)
{
  // ----------------------------------------------------------
  // Pack A
  //
  // q = r0 + i, counted from the bottom.
  // ----------------------------------------------------------

  for (int k = 0; k < srcWidth; ++k)
  {
    const int end = srcEnd[k];

    for (int i = 0; i < M; ++i)
    {
      const int q = r0 + i;

      A[i * srcWidth + k] =
        l[end - q];
    }
  }

  // ----------------------------------------------------------
  // Gather FC
  //
  // For target p = p0 + j, only q <= done + p
  // belongs to the current reverse update.
  // ----------------------------------------------------------

  for (int j = 0; j < N; ++j)
  {
    const int p = p0 + j;

    const int target =
      rowindices[eK - 1 - done - p];

    const int ref_pos =
      colpointers[target + 1] - 1;

    for (int i = 0; i < M; ++i)
    {
      const int q = r0 + i;

      if (q <= done + p)
      {
        const int ndx =
          rowindices[eK - 1 - q];

        const int pos =
          ref_pos - indmap[ndx];

        FC[i * N + j] =
          f[pos];
      }
      else
      {
        // Entry was not part of the forward C block.
        FC[i * N + j] = 0.0;
      }
    }
  }

  // ----------------------------------------------------------
  // Reverse dense multiplication
  //
  //     FA  = FC Bt^T
  //     FBt = A^T FC
  // ----------------------------------------------------------

  ADmatmul_block<M, N>(
      A.data(),
      Bt.data(),
      FC.data(),
      FA.data(),
      FBt.data(),
      srcWidth,  // lda
      N,         // ldbt
      N,         // ldfc
      srcWidth,  // ldfa
      N,         // ldfbt
      srcWidth   // K
  );

  // ----------------------------------------------------------
  // Scatter FA back to source supernode K
  // ----------------------------------------------------------

  for (int i = 0; i < M; ++i)
  {
    const int q = r0 + i;

    for (int k = 0; k < srcWidth; ++k)
    {
      const int pos =
        srcEnd[k] - q;

      f[pos] -=
        FA[i * srcWidth + k];
    }
  }

  // ----------------------------------------------------------
  // Scatter FBt back to source supernode K
  // ----------------------------------------------------------

  for (int j = 0; j < N; ++j)
  {
    const int q =
      done + p0 + j;

    for (int k = 0; k < srcWidth; ++k)
    {
      const int pos =
        srcEnd[k] - q;

      f[pos] -=
        FBt[k * N + j];
    }
  }
}

#endif
