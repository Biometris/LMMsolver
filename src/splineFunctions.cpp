#include <RcppArmadillo.h>
// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Row-wise kronecker product
//'
//' Row-wise kronecker product
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

//' Spectral decomposition of D'D, returns a q x (q-ord) matrix
//'
//' Spectral decomposition of D'D, returns a q x (q-ord) matrix
//'
//' @export
// [[Rcpp::export]]
arma::mat calcUsc(const double& q,
                  const double& ord) {
  arma::mat D = arma::diff(arma::eye(q, q), ord);
  arma::mat DtD = D.t() * D;

  arma::vec eigVals( D.n_cols );
  arma::mat eigVecs( size(D) );
  arma::eig_sym(eigVals, eigVecs, DtD);
  // fliplr and reverse are needed because eigenVals and eigenVecs are returned
  // sorted in ascending order.
  eigVecs = fliplr( eigVecs.tail_cols(q - ord) );
  arma::mat S = arma::diagmat( 1 / sqrt(reverse(eigVals.tail(q - ord))) );
  return eigVecs * S;
}

