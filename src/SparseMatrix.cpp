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
Rcpp::S4 RowKron(Rcpp::S4 sX1, Rcpp::S4 sX2)
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
  // see Masaki Tsuda, Rcpp for everyone,
  // https://teuder.github.io/rcpp4everyone_en/
  S4 L("spam");
  L.slot("entries") = entries;
  L.slot("colindices") = colindices;
  L.slot("rowpointers") = rowpointers;
  L.slot("dimension") = dimension;
  return L;
}

// help function, to count the number of non-zero's in C=A*B
int cntProduct(const SparseMatrix& A, const SparseMatrix B)
{
  const int nRows = A.dim[0];

  int cnt=0;
  vector<bool> S(B.dim[1],false);
  vector<int> ndx(B.dim[1],-1);
  for (int i=0;i<nRows;i++)
  {
    int ptr_ndx = 0;
    int s1 = A.rowpointers[i];
    int e1 = A.rowpointers[i+1];
    for (int j=s1;j<e1;j++)
    {
      int m = A.colindices[j];
      int s2 = B.rowpointers[m];
      int e2 = B.rowpointers[m+1];
      for (int k=s2;k<e2;k++)
      {
        int index = B.colindices[k];
        if(!S[index])
        {
          S[index] = true;
          cnt++;
          ndx[ptr_ndx++] = index;
        }
      }
    }
    for (int j=0;j<ptr_ndx;j++)
    {
      S[ndx[j]] = false;
      ndx[j] = -1;
    }
  }
  return cnt;
}


// [[Rcpp::export]]
Rcpp::S4 MatrixProduct(Rcpp::S4 sA, Rcpp::S4 sB)
{
  if ((as<std::string>(sA.attr("class")) != "spam") ||
      (as<std::string>(sB.attr("class")) != "spam"))
  {
    stop("Both arguments for MatrixProduct should be of class spam");
  }
  SparseMatrix A(sA);
  SparseMatrix B(sB);
  if (A.dim[1] != B.dim[0])
  {
    stop("MatrixProduct wrong dimensions");
  }
  const int nRows = A.dim[0];
  const int nCols = B.dim[1];
  const int N = cntProduct(A, B);

  IntegerVector dimension(2);
  IntegerVector rowpointers(nRows+1);
  rowpointers[nRows] = N+1;
  dimension[0] = nRows;
  dimension[1] = nCols;

  IntegerVector colindices(N,-1);
  NumericVector entries(N, 0.0);

  int cnt=0;
  vector<bool> S(B.dim[1],false);
  vector<double> val(B.dim[1],0.0);
  vector<int> ndx(B.dim[1],-1);
  for (int i=0;i<nRows;i++)
  {
    int ptr_ndx = 0;
    rowpointers[i] = cnt + 1;
    int s1 = A.rowpointers[i];
    int e1 = A.rowpointers[i+1];
    for (int j=s1;j<e1;j++)
    {
      double alpha = A.entries[j];
      int m = A.colindices[j];
      int s2 = B.rowpointers[m];
      int e2 = B.rowpointers[m+1];
      for (int k=s2;k<e2;k++)
      {
        int index = B.colindices[k];
        if (!S[index]) {
          S[index] = true;
          ndx[ptr_ndx++] = index;
        }
        val[index] += alpha*B.entries[k];
      }
    }
    sort(ndx.begin(),ndx.begin()+ptr_ndx);
    for (int j=0;j<ptr_ndx;j++) {
      colindices[cnt] = ndx[j] + 1;
      entries[cnt] = val[ndx[j]];
      cnt++;
      S[ndx[j]] = false;
      val[ndx[j]] = 0.0;
      ndx[j] = -1;
    }
  }

  // see Masaki Tsuda, Rcpp for everyone,
  // https://teuder.github.io/rcpp4everyone_en/
  S4 L("spam");
  L.slot("entries") = entries;
  L.slot("colindices") = colindices;
  L.slot("rowpointers") = rowpointers;
  L.slot("dimension") = dimension;
  return L;
}


/*
 * Not used, only for diagnostics
 *
 // [[Rcpp::export]]
 int cntProduct(Rcpp::S4 sA, const Rcpp::S4 sB)
 {
 SparseMatrix A(sA);
 SparseMatrix B(sB);
 int cnt = cntProduct(A, B);
 return cnt;
 }

 */

