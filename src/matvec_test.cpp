#include <Rcpp.h>
#include <vector>

#include <chrono>

#include "matvec_test.h"

using namespace std::chrono;

using namespace Rcpp;
using namespace std;


inline void matmulImplNaive(const double* left,
                            const double* right,
                            double* result,
                            int rows,
                            int columns,
                            int inners)
{
  for (int row = 0; row < rows; row++) {
    for (int col = 0; col < columns; col++) {
      for (int inner = 0; inner < inners; inner++) {
        result[row * columns + col] +=
          left[row * inners + inner] *
          right[inner * columns + col];
      }
    }
  }
}


inline void matmulImplLoopOrder(const double *left,
                                const double *right,
                                double *result,
                                int rows,
                                int columns,
                                int inners)
{
  for (int row = 0; row < rows; row++) {
    for (int inner = 0; inner < inners; inner++) {
      for (int col = 0; col < columns; col++) {
        result[row * columns + col] +=
          left[row * inners + inner] *
          right[inner * columns + col];
      }
    }
  }
}

inline void matmulImplTiling(const double *left, const double *right,
                             double *result, int rows, int columns, int inners, int tileSize) {
  for (int innerTile = 0; innerTile < inners; innerTile += tileSize) {
    for (int row = 0; row < rows; row++) {
      int innerTileEnd = std::min(inners, innerTile + tileSize);
      for (int inner = innerTile; inner < innerTileEnd; inner++) {
        for (int column = 0; column < columns; column++) {
          result[row * columns + column] +=
            left[row * inners + inner] * right[inner * columns + column];
        } } } } }



#include <Rcpp.h>
#include <vector>
#include <chrono>

using namespace Rcpp;
using namespace std::chrono;


// ------------------------------------------------------------
// 4 x 4 micro-kernel
//
// Computes:
//
//     C4 = A4 %*% B4
//
// where:
//   A4 = 4 x K
//   B4 = K x 4
//
// All matrices are stored in row-major order.
// ------------------------------------------------------------

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


// ------------------------------------------------------------
// Rcpp interface
//
// Computes:
//
//     C = A %*% B
//
// using 4 x 4 output blocks.
// ------------------------------------------------------------

// [[Rcpp::export]]
NumericMatrix matmul4x4_blocks(
    NumericMatrix A,
    NumericMatrix B,
    int reps)
{
  const int rows = A.nrow();
  const int K    = A.ncol();

  if (B.nrow() != K)
    stop("Incompatible dimensions: ncol(A) must equal nrow(B)");

  const int cols = B.ncol();

  if (rows % 4 != 0 || cols % 4 != 0)
    stop("Number of rows of A and columns of B must be divisible by 4");

  // ----------------------------------------------------------
  // Convert R column-major matrices to row-major storage.
  // ----------------------------------------------------------

  std::vector<double> left(rows * K);
  std::vector<double> right(K * cols);
  std::vector<double> result(rows * cols);

  for (int i = 0; i < rows; ++i)
  {
    for (int k = 0; k < K; ++k)
    {
      left[i * K + k] = A(i, k);
    }
  }

  for (int k = 0; k < K; ++k)
  {
    for (int j = 0; j < cols; ++j)
    {
      right[k * cols + j] = B(k, j);
    }
  }

  // ----------------------------------------------------------
  // Benchmark.
  // ----------------------------------------------------------

  auto start = high_resolution_clock::now();

  for (int r = 0; r < reps; ++r)
  {
    std::fill(result.begin(), result.end(), 0.0);

    for (int i = 0; i < rows; i += 4)
    {
      for (int j = 0; j < cols; j += 4)
      {
        matmul4x4(
          &left[i * K],
               &right[j],
               &result[i * cols + j],
               K,
               cols,
               cols,
               K
        );
      }
    }
  }

  auto end = high_resolution_clock::now();

  double elapsed =
    duration<double>(end - start).count() / reps;

  Rcpp::Rcout
  << "Average time: "
  << elapsed
  << " seconds\n";

  // ----------------------------------------------------------
  // Convert row-major result back to R column-major matrix.
  // ----------------------------------------------------------

  NumericMatrix out(rows, cols);

  for (int i = 0; i < rows; ++i)
  {
    for (int j = 0; j < cols; ++j)
    {
      out(i, j) = result[i * cols + j];
    }
  }

  return out;
}


// [[Rcpp::export]]
NumericMatrix matmul_test(NumericMatrix A, NumericMatrix B, int reps, int tileSize)
{
  const int rows    = A.nrow();
  const int inners  = A.ncol();
  const int columns = B.ncol();

  if (B.nrow() != inners)
    stop("Incompatible dimensions");

  std::vector<double> left(rows * inners);
  std::vector<double> right(inners * columns);
  std::vector<double> result(rows * columns, 0.0);

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < inners; j++)
      left[i * inners + j] = A(i, j);

  for (int i = 0; i < inners; i++)
    for (int j = 0; j < columns; j++)
      right[i * columns + j] = B(i, j);

  auto start = std::chrono::high_resolution_clock::now();

  for (int r = 0; r < reps; r++) {
    std::fill(result.begin(), result.end(), 0.0);

    if (tileSize == -1) {
      matmulImplNaive(
        left.data(), right.data(), result.data(),
        rows, columns, inners
      );
    } else if (tileSize == 0) {
      matmulImplLoopOrder(
        left.data(), right.data(), result.data(),
        rows, columns, inners
      );
    } else {
      matmulImplTiling(
        left.data(), right.data(), result.data(),
        rows, columns, inners, tileSize
      );
    }
  }

  auto end = std::chrono::high_resolution_clock::now();

  double elapsed = std::chrono::duration<double>(end - start).count() / reps;

  Rcpp::Rcout << "Average time: " << elapsed << " seconds\n";

  // Convert result vector to an R matrix
  NumericMatrix out(rows, columns);

  for (int i = 0; i < rows; i++)
    for (int j = 0; j < columns; j++)
      out(i, j) = result[i * columns + j];

  return out;
}

