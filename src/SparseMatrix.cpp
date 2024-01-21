#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include "AuxFun.h"
#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(Rcpp::S4 obj)
{
  if (as<std::string>(obj.attr("class")) != "spam") {
    std::string str = "wrong class " +
      as<std::string>(obj.attr("class")) +
      " for SparseMatrix constructor, should be class spam.";
    stop(str);
  }
  // get numeric values:
  entries = GetNumericVector(obj, "entries");
  // use 0-based indexing in C++/Rcpp:
  colindices = GetIntVector(obj, "colindices", 0);
  rowpointers = GetIntVector(obj, "rowpointers", 0);
  // dimensions of matrix:
  dim = GetIntVector(obj, "dimension", 1);
}

// [[Rcpp::export]]
List RowKron(Rcpp::S4 sX1, Rcpp::S4 sX2)
{
  SparseMatrix X1(sX1);
  SparseMatrix X2(sX2);
  const int n = X1.dim[0];
  const int q1 = X1.dim[1];
  const int q2 = X2.dim[1];
  const int q = q1*q2;
  IntegerVector dimension(2);
  IntegerVector rowpointers(n+1);
  dimension[0] = n;
  dimension[1] = q;
  int N = 0;
  for (int r=0;r<n;r++) {
    rowpointers[r] = N+1;
    int n1 = X1.rowpointers[r+1]-X1.rowpointers[r];
    int n2 = X2.rowpointers[r+1]-X2.rowpointers[r];
    N += n1*n2;
  }
  rowpointers[n] = N+1;

  IntegerVector colindices(N);
  NumericVector entries(N);

  int k = 0;
  for (int r=0;r<n;r++) {
    int s1 = X1.rowpointers[r];
    int e1 = X1.rowpointers[r+1];

    int s2 = X2.rowpointers[r];
    int e2 = X2.rowpointers[r+1];
    for (int j1=s1;j1<e1;j1++) {
      for (int j2=s2;j2<e2;j2++) {
        int ndx1 = X1.colindices[j1];
        int ndx2 = X2.colindices[j2];
        colindices[k] = ndx2 + q2*ndx1+1;
        entries[k] = X1.entries[j1]*X2.entries[j2];
        k++;
      }
    }
  }
  List L = List::create(Named("entries") = entries,
                        Named("colindices") = colindices,
                        Named("rowpointers") = rowpointers,
                        Named("dimension") = dimension);
  return L;
}

