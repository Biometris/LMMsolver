#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

#include "AuxFun.h"
#include "SparseMatrix.h"

SparseMatrix::SparseMatrix(Rcpp::S4 obj)
{
  entries = Rcpp::clone<Rcpp::NumericVector>(obj.slot("entries"));
  colindices = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colindices"));
  rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("rowpointers"));
  dim = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("dimension"));
  // in C/C++ indices start from zero, in R from one:
  transf2C(colindices);
  transf2C(rowpointers);
}

// [[Rcpp::export]]
List RowKron(SEXP sX1, SEXP sX2)
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

  //IntegerVector row(N);
  IntegerVector colindices(N);
  NumericVector entries(N);

  int k = 0;
  for (int r=0;r<n;r++) {
    //Rcout << "row " << r << endl;
    int s1 = X1.rowpointers[r];
    int e1 = X1.rowpointers[r+1];

    int s2 = X2.rowpointers[r];
    int e2 = X2.rowpointers[r+1];
    for (int j1=s1;j1<e1;j1++) {
      for (int j2=s2;j2<e2;j2++) {
        int ndx1 = X1.colindices[j1];
        int ndx2 = X2.colindices[j2];
        //Rcout << "ndx " << ndx2 + q2*ndx1 << " val "
        //      << X1.entries[j1]*X2.entries[j2] << endl;
        //row[k] = r+1;
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

/*

Not used at the moment: functions to calculate GLAM:

NumericVector prod(const SparseMatrix& M, const NumericVector& y, int start_y)
{
  NumericVector result(M.dim[0]);
  for (int r=0;r<M.dim[0];r++)
  {
    int s = M.rowpointers[r];
    int e = M.rowpointers[r+1];
    double z = 0.0;
    for (int k=s;k<e;k++)
    {
      double v = M.entries[k];
      int col = M.colindices[k];
      z += v*y[start_y+col];
    }
    result[r] = z;
  }
  return result;
}

NumericVector prodAux(const SparseMatrix& B, int nBlk, NumericVector x)
{
  int n = B.dim[0];
  int q = B.dim[1];

  NumericVector z(nBlk*n);
  int s = 0;
  for (int i=0;i<nBlk;i++)
  {
    NumericVector tmp = prod(B, x, s);
    s += q;
    for (int j=0;j<n;j++)
    {
      z[nBlk*j+i] = tmp[j];
    }
  }
  return z;
}

// [[Rcpp::export]]
NumericVector KronProd2(SEXP C1, SEXP C2, const NumericVector& y)
{
  SparseMatrix B1(C1);
  SparseMatrix B2(C2);
  int q1 = B1.dim[1];
  int n2 = B2.dim[0];

  NumericVector z1 = prodAux(B2, q1, y);
  NumericVector z2 = prodAux(B1, n2, z1);
  return z2;
}

// [[Rcpp::export]]
NumericVector KronProd(SEXP C1, SEXP C2, SEXP C3, const NumericVector& y)
{
  SparseMatrix B1(C1);
  SparseMatrix B2(C2);
  SparseMatrix B3(C3);
  int q1 = B1.dim[1];
  int q2 = B2.dim[1];
  int n2 = B2.dim[0];
  int n3 = B3.dim[0];

  NumericVector z1 = prodAux(B3, q1*q2, y);
  NumericVector z2 = prodAux(B2, q1*n3, z1);
  NumericVector z3 = prodAux(B1, n2*n3, z2);
  return z3;
}


// [[Rcpp::export]]
NumericVector KronProdList(List L, const NumericVector& y)
{
  int sz = L.size();
  vector<SparseMatrix> B;
  for (int i=0;i<sz;i++)
  {
    Rcpp::S4 obj = L[i];
    B.push_back(obj);
  }
  NumericVector z = y;
  for (int i=sz-1;i>=0;i--)
  {
    int nBlk = 1;
    for (int k=0;k<sz;k++)
    {
      if (k < i) nBlk *= B[k].dim[1];
      if (k > i) nBlk *= B[k].dim[0];
    }
    z = prodAux(B[i], nBlk, z);
  }
  return z;
}

// [[Rcpp::export]]
IntegerVector getOrder(SEXP A, int q1, int q2,
                        const IntegerVector& s1,
                        const IntegerVector& s2)
{
  SparseMatrix M(A);
  int nElem = M.entries.size();
  IntegerVector x(nElem);
  int dim2 = s2.size();

  map<int,int> Key1, Key2;
  for (int i=0;i<s1.size();i++) {
    Key1[s1[i]-1] = i;
  }
  for (int i=0;i<s2.size();i++) {
    Key2[s2[i]-1] = i;
  }

  for (int i=0;i<M.dim[0];i++)
  {
    int s = M.rowpointers[i];
    int e = M.rowpointers[i+1];

    int i1 = i / q2;
    int i2 = i % q2;

    for (int k=s;k<e;k++)
    {
      int j = M.colindices[k];
      int j1 = j / q2;
      int j2 = j % q2;

      int ndx1 = j1*q1 + i1;
      int ndx2 = j2*q2 + i2;
      int k1 = Key1[ndx1];
      int k2 = Key2[ndx2];

      // plus one because of the R-index:
      x[k] = k1*dim2 + k2 + 1;
    }
  }
  return x;
}

*/

