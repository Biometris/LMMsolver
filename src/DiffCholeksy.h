#ifndef CALC_DIFFERENTIATION_CHOLESKY_ALGORITHM
#define CALC_DIFFERENTIATION_CHOLESKY_ALGORITHM

#include <Rcpp.h>

using namespace Rcpp;
using namespace std;

// auxiliary class for function logdet
class SparseMatrix
{
public:
  SparseMatrix(const IntegerVector& colpointers,
               const IntegerVector& rowindices,
               const NumericVector& entries) : _colptr(colpointers),
                                               _rowndx(rowindices),
                                               _entries(entries)
  { }
  ~SparseMatrix() {}
  int getpos(int r, int c) const
  {
    int s = _colptr[c];
    int e = _colptr[c+1];
    for (int k=s;k<e;k++)
    {
      if (_rowndx[k] == r) return k;
    }
    return -1;
  }
  double getvalue(int r, int c) const
  {
    int p = getpos(r,c);
    if (p==-1) return 0.0;
    return _entries[p];
  }
  double&  operator[](int k) { return _entries[k]; }
  const double& operator[](int k) const { return _entries[k]; }
  int rowndx(int k) const { return _rowndx[k]; }

  pair<int,int> getcol(int c)
  {
    int s = _colptr[c];
    int e = _colptr[c+1];
    return pair<int,int>(s,e);
  }
private:
  const IntegerVector _colptr;
  const IntegerVector _rowndx;
  NumericVector _entries;
};

#endif

