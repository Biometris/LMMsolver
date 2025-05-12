// sparse left-looking Cholesky using supernodes.
//
// For details see:
// Ng, Esmond G., and Barry W. Peyton.,
// "Block sparse Cholesky algorithms on advanced uniprocessor computers."
// SIAM Journal on Scientific Computing 14, no. 5 (1993): 1034-1056.
//
// Furrer, Reinhard, and Stephan R. Sain.
// "spam: A sparse matrix R package with emphasis on MCMC
// methods for Gaussian Markov random fields."
// Journal of Statistical Software 36 (2010): 1-25.

#include <Rcpp.h>
#include <set>
#include <vector>
#include "AuxFun.h"
#include "SparseMatrix.h"
#include "cholesky.h"

using namespace Rcpp;
using namespace std;

// make indmap for supernode J:
void makeIndMap(IntegerVector& indmap,
                int J,
                const IntegerVector& rowpointers,
                const IntegerVector& rowindices)
{
  int s = rowpointers[J];
  int e = rowpointers[J+1];
  int l = 0;
  for (int i=e-1; i>=s;i--)
  {
    indmap[rowindices[i]] = l++;
  }
}

// j is current column in Supernode J
void cmod1(NumericVector& L, int j, int J,
           const IntegerVector& supernodes,
           const IntegerVector& colpointers)
{
  const int& s = colpointers[j];
  const int& e = colpointers[j+1];
  // for all columns in supernode J left to j:
  for (int k=supernodes[J];k<j;k++)
  {
    const int& jk = colpointers[k] + (j-k);
    int ik = jk;
    const double& Ljk = L[jk];
    for (int ij=s; ij<e; ij++)
    {
       L[ij] -= L[ik++]*Ljk;
       //ik++;
    }
  }
}

// Adjust column j for all columns in supernode K:
void cmod2(NumericVector& L, int j, int K, int sz,
           NumericVector& t,
           const IntegerVector& indmap,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  // init t:
  for (int i=0;i<sz;i++)
  {
    t[i] = 0.0;
  }

  const int& sCol = supernodes[K];
  const int& eCol = supernodes[K+1];

  for (int k=sCol; k<eCol; k++)
  {
    int jk = colpointers[k+1]-sz;
    int ik = jk;
    const double& Ljk = L[jk];
    for (int i=sz-1;i>=0;i--)
    {
      t[i] += L[ik++]*Ljk;
      //ik++;
    }
  }

  int r = rowpointers[K+1]-1;
  int ref_pos = colpointers[j+1] - 1;
  for (int i=0;i<sz;i++)
  {
    int ndx = rowindices[r--];
    int pos = ref_pos - indmap[ndx];
    L[pos] -= t[i];
  }
}

void cdiv(NumericVector& L, int j, const IntegerVector& colpointers)
{
  const int& s = colpointers[j];
  const int& e = colpointers[j+1];

  // pivot:
  L[s] = sqrt(L[s]);
  // update column j:
  double Ls = L[s];
  for (int i = s + 1; i < e; i++)
  {
    L[i] /= Ls;
  }
}


void cholesky(NumericVector& L,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  const int N = colpointers.size() - 1;
  const int Nsupernodes = supernodes.size()-1;

  // linked lists, see section 4.2 Ng and Peyton
  IntegerVector HEAD(N,-1);
  IntegerVector LINK(Nsupernodes,-1);

  IntegerVector colhead = clone(rowpointers);
  for (int J=0; J<Nsupernodes;J++)
  {
    int szNode = supernodes[J+1] - supernodes[J];
    colhead[J] += szNode-1;
    if (colhead[J] < rowpointers[J+1]-1)
    {
      int rNdx = rowindices[colhead[J]+1];
      insert(HEAD, LINK, rNdx, J);
    }
  }

  IntegerVector indmap(N,0);
  NumericVector t(N);
  for (int J=0; J<Nsupernodes;J++) {
    makeIndMap(indmap, J, rowpointers, rowindices);
    for (int j=supernodes[J];j<supernodes[J+1];j++)
    {
      int K = HEAD[j];
      while (K!=-1)
      {
        int nextK = LINK[K];
        int sz = rowpointers[K+1] - colhead[K];
        cmod2(L, j, K, sz, t, indmap, supernodes, rowpointers, colpointers, rowindices);

        colhead[K]++;
        if (colhead[K] < rowpointers[K+1])
        {
          int rNdx = rowindices[colhead[K]];
          insert(HEAD, LINK, rNdx, K);
        }
        K = nextK;
      }
      HEAD[j] = -1;
      cmod1(L, j, J, supernodes, colpointers);
      cdiv(L, j, colpointers);
    }
    colhead[J]++;
  }
}

double logdet(const NumericVector& L, const IntegerVector& colpointers)
{
  const int N = colpointers.size() - 1;
  double sum = 0;
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    sum += 2.0*log(L[s]);
  }
  return sum;
}


NumericVector forwardCholesky(
    const NumericVector& L,
    const NumericVector& b,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices,
    const IntegerVector& pivot,
    const IntegerVector& invpivot)
{
  const int Nsupernodes = supernodes.size()-1;
  const int N = colpointers.size() - 1;
  NumericVector x(N);
  NumericVector Pb(N);  // permutation of P
  NumericVector sum(N); // sum for each column
  for (int i=0;i<N;i++)
  {
    Pb[i] = b[pivot[i]];
  }
  for (int J=0; J<Nsupernodes;J++)
  {
    int s = rowpointers[J];
    for (int j=supernodes[J]; j<supernodes[J+1]; j++)
    {
      const double x_j = (Pb[j]-sum[j])/L[colpointers[j]];
      x[j] = x_j;
      // the non-diagonal elements of column j
      for (int ndx = colpointers[j]+1, k=s+1; ndx < colpointers[j+1]; ndx++)
      {
        int i = rowindices[k++];
        sum[i] += L[ndx]*x_j;
      }
      s++;
    }
  }

  NumericVector xP(N); // inverse permutation
  for (int i=0;i<N;i++) {
    xP[i] = x[invpivot[i]];
  }
  return xP;

}


NumericVector backwardCholesky(
    const NumericVector& L,
    const NumericVector& b,
    const IntegerVector& supernodes,
    const IntegerVector& rowpointers,
    const IntegerVector& colpointers,
    const IntegerVector& rowindices,
    const IntegerVector& pivot,
    const IntegerVector& invpivot)
{
  const int Nsupernodes = supernodes.size()-1;
  const int N = colpointers.size() - 1;
  NumericVector x(N);
  NumericVector Pb(N);  // permutation of P
  NumericVector sum(N); // sum for each column
  for (int i=0;i<N;i++)
  {
    Pb[i] = b[pivot[i]];
  }
  for (int J=Nsupernodes-1; J>=0;J--)
  {
    int NodeSz = supernodes[J+1] - supernodes[J];
    int s = rowpointers[J]+(NodeSz-1);
    for (int j=supernodes[J+1]-1; j>=supernodes[J]; j--)
    {
      //Rcout << "column " << j << endl;
      double alpha = L[colpointers[j]];
      double x_j = Pb[j];
      // the non-diagonal elements of column j
      for (int ndx = colpointers[j]+1, k=s+1; ndx < colpointers[j+1]; ndx++)
      {
        int i = rowindices[k++];
        x_j -= L[ndx]*x[i];
      }
      x[j] = x_j/alpha;
      s--;
    }
  }

  NumericVector xP(N); // inverse permutation
  for (int i=0;i<N;i++) {
    xP[i] = x[invpivot[i]];
  }
  return xP;
}


/*
 *
// [[Rcpp::export]]
double logdet(Rcpp::S4 obj, NumericVector lambda)
{
  //Rcpp::S4 obj(arg);
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericVector L = obj.slot("entries");
  NumericMatrix P = obj.slot("P");
  //NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

  // define matrix L (lower triangle matrix values)
  const int sz = P.nrow();
  const int n_prec_mat = P.ncol();
  //NumericVector L(sz, 0.0);
  for (int i=0;i<sz;i++)
  {
    L[i] = 0.0;
  }
  for (int k=0;k<n_prec_mat;k++)
  {
    NumericMatrix::Column Pk = P(_, k);
    double alpha = lambda[k];
    for (int i=0;i<sz;i++)
    {
      L[i] += alpha*Pk[i];
    }
  }
  cholesky(L, supernodes, rowpointers, colpointers, rowindices);
  return logdet(L, colpointers);
}

*/

/*
// Just to show the structure of the sparse cholesky matrix with supernodes.
// [[Rcpp::export]]
NumericMatrix PrintCholesky(Rcpp::S4 obj)
{
  Rcout << "Class: " << as<std::string>(obj.attr("class")) << std::endl;

  IntegerVector supernodes = GetIntVector(obj, "supernodes", 0);
  // Exchange row and columns compared to spam object, as in Ng and Peyton 1993
  IntegerVector colpointers = GetIntVector(obj, "rowpointers", 0);
  IntegerVector rowpointers = GetIntVector(obj, "colpointers", 0);
  IntegerVector rowindices = GetIntVector(obj, "colindices", 0);
  IntegerVector pivot = GetIntVector(obj, "pivot", 0);
  IntegerVector invpivot = GetIntVector(obj, "invpivot", 0);

  NumericVector L = Rcpp::clone<Rcpp::NumericVector>(obj.slot("entries"));

  const int Nsupernodes = supernodes.size()-1;
  const int N = colpointers.size() - 1;
  NumericMatrix A(N, N);
  for (int J=0; J<Nsupernodes;J++)
  {
    int s = rowpointers[J];
    Rcout << "Supernode: " << J << endl;
    for (int j=supernodes[J]; j<supernodes[J+1]; j++)
    {
      Rcout << "  Column: " << j << endl;
      int k = s;
      for (int ndx = colpointers[j]; ndx < colpointers[j+1]; ndx++)
      {
        int i = rowindices[k++];
        Rcout << "    row: " << i << " (ndx or key " << ndx << ")" << endl;
        A(i, j) = L[ndx];
      }
      s++;
    }
  }
  return A;
}

*/
