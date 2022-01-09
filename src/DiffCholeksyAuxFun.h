#ifndef DIFFCHOLESKYAUXFUN_HEADER
#define DIFFCHOLESKYAUXFUN_HEADER


#include <Rcpp.h>
#include <set>
#include <vector>

using namespace Rcpp;
using namespace std;

// Transform to C++ Notation indices
void transf2C(IntegerVector& ndx);

// not very efficient (but not too bad): Make a Class?
double getvalueC(IntegerVector rowpointers,
                 IntegerVector colindices,
                 NumericVector entries,
                 int i,
                 int j);

vector<double> convert_matrix(const vector<int>& rowindices_ext,
                              const vector<int>& colindices_ext,
                              const IntegerVector& pivot,
                              SEXP A);

#endif
