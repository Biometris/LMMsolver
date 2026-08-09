// Backwards Automatic Differentiation of the Cholesky Algorithm
// for calculating partial derivatives of the log-determinant of
// positive definite symmetric sparse matrices.
//
// Implementation by Martin Boer, 2026.
//
// The algorithm combines the sparse Cholesky factorization of
// Ng and Peyton (1993) with the reverse differentiation approach
// of Smith (1995).
//
// References:
//
// Ng, E. G. and B. W. Peyton (1993).
// "Block sparse Cholesky algorithms on advanced uniprocessor computers."
// SIAM Journal on Scientific Computing 14, 1034-1056.
//
// Furrer, R. and S. R. Sain (2010).
// "spam: A sparse matrix R package with emphasis on MCMC
// methods for Gaussian Markov random fields."
// Journal of Statistical Software 36, 1-25.
//
// Smith, S. P. (1995).
// "Differentiation of the Cholesky algorithm."
// Journal of Computational and Graphical Statistics 4, 134-147.
//
// Smith, S. P. (2000).
// "A Tutorial on Simplicity and Computational Differentiation
// for Statisticians."

#ifndef ADCHOLESKY_HEADER
#define ADCHOLESKY_HEADER

#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

void ADcmod1(NumericVector& F,
             const NumericVector& L, int j, int J,
             const IntegerVector& supernodes,
             const IntegerVector& colpointers);

void ADcmod2(NumericVector& F,
             const NumericVector& L, int j, int K, int sz,
             NumericVector& t,
             const IntegerVector& indmap,
             const IntegerVector& supernodes,
             const IntegerVector& rowpointers,
             const IntegerVector& colpointers,
             const IntegerVector& rowindices);

void ADcdiv(NumericVector& F,
            const NumericVector& L, int j, const IntegerVector& colpointers);

void ADcholesky(NumericVector& F,
                  const NumericVector& L,
                  const IntegerVector& supernodes,
                  const IntegerVector& rowpointers,
                  const IntegerVector& colpointers,
                  const IntegerVector& rowindices);

void initAD(NumericVector& F, const NumericVector& L, const IntegerVector& colpointers);


#endif
