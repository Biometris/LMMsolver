#ifndef DIFFCHOLESKYAUXFUN_HEADER
#define DIFFCHOLESKYAUXFUN_HEADER

#include <Rcpp.h>
#include <set>
#include <vector>

using namespace Rcpp;
using namespace std;

// Transform to C++ Notation indices
void transf2C(IntegerVector& ndx);

IntegerVector GetIntVector(Rcpp::S4 obj, const String& slotName, int ArrayIndexing);

NumericVector GetNumericVector(Rcpp::S4 obj, const String& slotName);

void insert(IntegerVector& HEAD, IntegerVector& LINK, int i, int J);

#endif
