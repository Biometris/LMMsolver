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

SparseMatrix::SparseMatrix(Rcpp::S4 obj)
{
  entries = Rcpp::clone<Rcpp::NumericVector>(obj.slot("entries"));
  colindices = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colindices"));
  rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("rowpointers"));
  dim = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("dimension"));
  // in C/C++ indices start from zero, in R from one:
  for (int i=0;i<colindices.size();i++) {
    colindices[i]--;
  }
  for (int i=0;i<rowpointers.size();i++) {
    rowpointers[i]--;
  }
}

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

//ReArrange <- function(A, q1, q2, s1, s2, z) {

// [[Rcpp::export]]
NumericVector ReArrange(SEXP A, int q1, int q2,
                        const IntegerVector& s1,
                        const IntegerVector& s2,
                        const NumericVector& z)
{
  SparseMatrix M(A);
  NumericVector x(z.size());
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
      x[k] = z[k1*dim2 + k2];
    }
  }
  return x;
}

