#include <Rcpp.h>
#include <set>
#include <vector>
#include "AuxFun.h"

using namespace Rcpp;
using namespace std;

// Transform to C++ Notation indices
void transf2C(IntegerVector& ndx)
{
  for (int i=0;i<ndx.size();i++)
  {
    ndx[i]--;
  }
}

// [[Rcpp::export]]
IntegerVector GetIntVector(Rcpp::S4 obj, const String& slotName, int ArrayIndexing)
{
  IntegerVector x = Rcpp::clone<Rcpp::IntegerVector>(obj.slot(slotName));
  if (ArrayIndexing == 0) {
    transf2C(x);
    return x;
  }
  if (ArrayIndexing == 1) {
    return x;
  }
  stop("argument ArrayIndex should be 0-based (C/C++) or 1-based (R).");
  return x;
}

NumericVector GetNumericVector(Rcpp::S4 obj, const String& slotName) {
  NumericVector x = Rcpp::clone<Rcpp::NumericVector>(obj.slot(slotName));
  return x;
}

