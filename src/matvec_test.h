#ifndef MATVEC_TEST_HEADER
#define MATVEC_TEST_HEADER

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

NumericMatrix matmul_test(NumericMatrix A, NumericMatrix B, int reps, int tileSize);

#endif
