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

