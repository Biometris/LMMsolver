#ifndef CHOLESKY_HEADER
#define CHOLESKY_HEADER

#include <Rcpp.h>
#include <set>
#include <vector>

using namespace Rcpp;
using namespace std;

void makeIndMap(IntegerVector& indmap,
                int J,
                const IntegerVector& rowpointers,
                const IntegerVector& rowindices);

void cholesky(NumericVector& L,
                const IntegerVector& supernodes,
                const IntegerVector& rowpointers,
                const IntegerVector& colpointers,
                const IntegerVector& rowindices);

double logdet(const NumericVector& L, const IntegerVector& colpointers);

NumericVector backwardCholesky(
    const NumericVector& L,
    const NumericVector& b,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices,
    const IntegerVector& pivot,
    const IntegerVector& invpivot);

NumericVector forwardCholesky(
    const NumericVector& L,
    const NumericVector& b,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices,
    const IntegerVector& pivot,
    const IntegerVector& invpivot);


#endif
