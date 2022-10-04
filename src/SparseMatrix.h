#ifndef SPARSEMATRIX_HEADER
#define SPARSEMATRIX_HEADER


#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class SparseMatrix {
public:
  SparseMatrix(Rcpp::S4 obj);
  NumericVector entries;
  IntegerVector colindices;
  IntegerVector rowpointers;
  IntegerVector dim;
};

#endif
