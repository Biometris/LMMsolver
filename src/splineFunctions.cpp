#include <RcppArmadillo.h>
// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Row-wise kronecker product
//'
//' Row-wise kronecker product
//'
//' @param X1 A matrix.
//' @param X2 A matrix.
//'
//' @return The row-wise kronecker product of X1 and X2.
//'
//' @export
// [[Rcpp::export]]
arma::mat RowKronecker(const arma::mat& X1,
                       const arma::mat& X2) {
  arma::mat one1 = arma::ones(1, X1.n_cols);
  arma::mat one2 = arma::ones(1, X2.n_cols);
  arma::mat rowKron = arma::kron(X1, one2) % arma::kron(one1, X2);
  return rowKron;
}



