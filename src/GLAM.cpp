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
  for (int i=0;i<colindices.size();i++) {
    colindices[i] -= 1;
  }
  for (int i=0;i<rowpointers.size();i++) {
    rowpointers[i] -= 1;
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

// [[Rcpp::export]]
NumericVector prod(SEXP C, NumericVector y)
{
  SparseMatrix M(C);
  NumericVector result = prod(M, y, 0);
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
NumericVector prod2(SEXP C1, SEXP C2, NumericVector y)
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
NumericVector prod3(SEXP C1, SEXP C2, SEXP C3, NumericVector y)
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
NumericVector prodList(List L, NumericVector y)
{
  int d = L.size();
  vector<SparseMatrix> B;
  for (int i=0;i<d;i++)
  {
    Rcpp::S4 obj = L[i];
    B.push_back(obj);
  }
  NumericVector z = y;
  for (int i=d-1;i>=0;i--)
  {
    int nBlk = 1;
    for (int k=0;k<d;k++)
    {
      if (k < i) nBlk *= B[k].dim[1];
      if (k > i) nBlk *= B[k].dim[0];
    }
    z = prodAux(B[i], nBlk, z);
  }
  return z;
}
