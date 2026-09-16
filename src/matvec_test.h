#ifndef MATVEC_TEST_HEADER
#define MATVEC_TEST_HEADER

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

NumericMatrix matmul_test(NumericMatrix A, NumericMatrix B, int reps, int tileSize);

inline void matmul4x4(
    const double* A,
    const double* B,
    double* C,
    int lda,
    int ldb,
    int ldc,
    int K)
{
  double c00 = 0.0, c01 = 0.0, c02 = 0.0, c03 = 0.0;
  double c10 = 0.0, c11 = 0.0, c12 = 0.0, c13 = 0.0;
  double c20 = 0.0, c21 = 0.0, c22 = 0.0, c23 = 0.0;
  double c30 = 0.0, c31 = 0.0, c32 = 0.0, c33 = 0.0;

  for (int k = 0; k < K; ++k)
  {
    // Four values from the current row of A
    const double a0 = A[0 * lda + k];
    const double a1 = A[1 * lda + k];
    const double a2 = A[2 * lda + k];
    const double a3 = A[3 * lda + k];

    // Four values from the current row of B
    const double b0 = B[k * ldb + 0];
    const double b1 = B[k * ldb + 1];
    const double b2 = B[k * ldb + 2];
    const double b3 = B[k * ldb + 3];

    c00 += a0 * b0;
    c01 += a0 * b1;
    c02 += a0 * b2;
    c03 += a0 * b3;

    c10 += a1 * b0;
    c11 += a1 * b1;
    c12 += a1 * b2;
    c13 += a1 * b3;

    c20 += a2 * b0;
    c21 += a2 * b1;
    c22 += a2 * b2;
    c23 += a2 * b3;

    c30 += a3 * b0;
    c31 += a3 * b1;
    c32 += a3 * b2;
    c33 += a3 * b3;
  }

  C[0 * ldc + 0] = c00;
  C[0 * ldc + 1] = c01;
  C[0 * ldc + 2] = c02;
  C[0 * ldc + 3] = c03;

  C[1 * ldc + 0] = c10;
  C[1 * ldc + 1] = c11;
  C[1 * ldc + 2] = c12;
  C[1 * ldc + 3] = c13;

  C[2 * ldc + 0] = c20;
  C[2 * ldc + 1] = c21;
  C[2 * ldc + 2] = c22;
  C[2 * ldc + 3] = c23;

  C[3 * ldc + 0] = c30;
  C[3 * ldc + 1] = c31;
  C[3 * ldc + 2] = c32;
  C[3 * ldc + 3] = c33;
}


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
// B x B matrix multiplication
//
// A  : B x K
// Bm : K x B
// C  : B x B
//
// C = A Bm
//
// B is a compile-time constant.
// ============================================================

template <int B>
inline void matmul_block_old(
    const double* A,
    const double* Bm,
    double* C,
    int lda,
    int ldb,
    int ldc,
    int K)
{
  double c[B][B] = {};

  for (int k = 0; k < K; ++k)
  {
    const double* b = Bm + k * ldb;

    #pragma GCC unroll B
    for (int i = 0; i < B; ++i)
    {
      const double a = A[i * lda + k];

      #pragma GCC unroll B
      for (int j = 0; j < B; ++j)
        c[i][j] += a * b[j];
    }
  }

  for (int i = 0; i < B; ++i)
  {
    for (int j = 0; j < B; ++j)
      C[i * ldc + j] = c[i][j];
  }
}

#endif
