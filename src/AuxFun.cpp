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
    ndx[i] -= 1;
  }
}

// not very efficient (but not too bad): Make a Class?
double getvalueC(IntegerVector rowpointers,
                 IntegerVector colindices,
                 NumericVector entries,
                 int i,
                 int j)
{
  int s = rowpointers[i];
  int e = rowpointers[i+1];
  for (int k=s;k<e;k++)
  {
    if (colindices[k] == j) return entries[k];
  }
  return 0.0;
}

vector<double> convert_matrix(const vector<int>& rowindices_ext,
                              const vector<int>& colindices_ext,
                              const IntegerVector& pivot,
                              SEXP A)
{
  Rcpp::S4 obj(A);
  IntegerVector rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("rowpointers"));
  IntegerVector colindices  = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colindices"));
  NumericVector entries     = obj.slot("entries") ;

  transf2C(rowpointers);
  transf2C(colindices);

  const int Nrow = rowindices_ext.size();
  std::vector<double> z(Nrow);
  for (int i=0;i<Nrow;i++)
  {
    int r = rowindices_ext[i];
    int c = colindices_ext[i];
    double x = getvalueC(rowpointers,colindices,entries,pivot[r],pivot[c]);
    z[i] = x;
  }
  return z;
}

