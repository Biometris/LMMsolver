#include <Rcpp.h>
#include <vector>
#include <chrono>

#include "matvec_test.h"

using namespace Rcpp;
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
NumericMatrix matmul4x4_blocks_template(
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
        matmul_block<4,4>(
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
NumericMatrix matmul2x2_blocks_template(
    NumericMatrix A,
    NumericMatrix B,
    int reps)
{
  const int rows = A.nrow();
  const int K    = A.ncol();

  if (B.nrow() != K)
    stop("Incompatible dimensions: ncol(A) must equal nrow(B)");

  const int cols = B.ncol();

  if (rows % 2 != 0 || cols % 2 != 0)
    stop("Number of rows of A and columns of B must be divisible by 2");

  // ----------------------------------------------------------
  // Convert R column-major matrices to row-major storage.
  // ----------------------------------------------------------

  std::vector<double> left(rows * K);
  std::vector<double> right(K * cols);
  std::vector<double> result(rows * cols);

  for (int i = 0; i < rows; ++i)
  {
    for (int k = 0; k < K; ++k)
      left[i * K + k] = A(i, k);
  }

  for (int k = 0; k < K; ++k)
  {
    for (int j = 0; j < cols; ++j)
      right[k * cols + j] = B(k, j);
  }

  // ----------------------------------------------------------
  // Benchmark.
  // ----------------------------------------------------------

  auto start = high_resolution_clock::now();

  for (int r = 0; r < reps; ++r)
  {
    std::fill(result.begin(), result.end(), 0.0);

    for (int i = 0; i < rows; i += 2)
    {
      for (int j = 0; j < cols; j += 2)
      {
        matmul_block<2,2>(
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
      out(i, j) = result[i * cols + j];
  }

  return out;
}


// [[Rcpp::export]]
NumericMatrix matmul8x8_blocks_template(
    NumericMatrix A,
    NumericMatrix B,
    int reps)
{
  const int rows = A.nrow();
  const int K    = A.ncol();

  if (B.nrow() != K)
    stop("Incompatible dimensions: ncol(A) must equal nrow(B)");

  const int cols = B.ncol();

  if (rows % 8 != 0 || cols % 8 != 0)
    stop("Number of rows of A and columns of B must be divisible by 8");

  // ----------------------------------------------------------
  // Convert R column-major matrices to row-major storage.
  // ----------------------------------------------------------

  std::vector<double> left(rows * K);
  std::vector<double> right(K * cols);
  std::vector<double> result(rows * cols);

  for (int i = 0; i < rows; ++i)
  {
    for (int k = 0; k < K; ++k)
      left[i * K + k] = A(i, k);
  }

  for (int k = 0; k < K; ++k)
  {
    for (int j = 0; j < cols; ++j)
      right[k * cols + j] = B(k, j);
  }

  // ----------------------------------------------------------
  // Benchmark.
  // ----------------------------------------------------------

  auto start = high_resolution_clock::now();

  for (int r = 0; r < reps; ++r)
  {
    std::fill(result.begin(), result.end(), 0.0);

    for (int i = 0; i < rows; i += 8)
    {
      for (int j = 0; j < cols; j += 8)
      {
        matmul_block<8,8>(
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
      out(i, j) = result[i * cols + j];
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

