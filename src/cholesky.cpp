// Backwards Automated Differentation of Cholesky Algorithm
// to calculate the partial derivatives of log-determinant of
// a positive definite symmetric (sparse) matrix.
//
// For details on the implementation of sparse Cholesky, see:
// Ng and Peyton 1993, Furrer and Sain 2010
//
// For details on Backwards Automated Differentiation see:
// S.P. Smith 1995, Differentiation of the Cholesky Algorithm
// and
// S.P. Smith 2000, A TUTORIAL ON SIMPLICITY AND COMPUTATIONAL DIFFERENTIATION FOR
// STATISTICIANS
//
#include <Rcpp.h>
#include <set>
#include <vector>
#include "AuxFun.h"
#include "NodeList.h"
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
  const int& sCol = supernodes[K];
  const int& eCol = supernodes[K+1];
  if (eCol - sCol > 1)
  {
    // init t:
    for (int i=0;i<sz;i++)
    {
      t[i] = 0.0;
    }

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
  } else {
    // if only one column in supernode:
    int k = sCol;
    int jk = colpointers[k+1]-sz;
    int ik = jk;
    const double& Ljk = L[jk];
    int r = rowpointers[K+1] - sz;
    int ref_pos = colpointers[j+1] - 1;
    for (int i=0;i<sz;i++)
    {
      int ndx = rowindices[r++];
      int pos = ref_pos - indmap[ndx];
      L[pos] -= L[ik++]*Ljk;
      //ik++;
    }
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
  vector<Node *> S(N);

  IntegerVector colhead = clone(rowpointers);
  for (int J=0; J<Nsupernodes;J++)
  {
    int szNode = supernodes[J+1] - supernodes[J];
    colhead[J] += szNode-1;
    if (colhead[J] < rowpointers[J+1]-1)
    {
      int rNdx = rowindices[colhead[J]+1];
      Node *nodeJ = new Node(J);
      S[rNdx] = add(nodeJ, S[rNdx]);
    }
  }
  IntegerVector indmap(N,0);
  NumericVector t(N);
  for (int J=0; J<Nsupernodes;J++) {
    makeIndMap(indmap, J, rowpointers, rowindices);
    for (int j=supernodes[J];j<supernodes[J+1];j++)
    {
      for (Node **ptr = &S[j]; *ptr;)
      {
        Node *f = removefirstnode(ptr);
        int K = f->data;
        int sz = rowpointers[K+1] - colhead[K];
        cmod2(L, j, K, sz, t, indmap, supernodes, rowpointers, colpointers, rowindices);

        colhead[K]++;
        if (colhead[K] < rowpointers[K+1]) {
          int rNdx = rowindices[colhead[K]];
          S[rNdx] = add(f, S[rNdx]);
        } else {
          delete f;
        }
      }
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

// Just to show the structure of the sparse cholesky matrix with supernodes.
// [[Rcpp::export]]
NumericMatrix PrintCholesky(SEXP cholC)
{
  Rcpp::S4 obj(cholC);
  // We use transpose for calculating Automated Differentiation.
  IntegerVector supernodes = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("supernodes"));
  IntegerVector colpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("rowpointers"));
  IntegerVector rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colpointers"));
  IntegerVector rowindices = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colindices"));
  IntegerVector pivot = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("pivot"));
  IntegerVector invpivot = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("invpivot"));

  NumericVector L = Rcpp::clone<Rcpp::NumericVector>(obj.slot("entries"));

  // C using indices starting at 0:
  transf2C(supernodes);
  transf2C(colpointers);
  transf2C(rowpointers);
  transf2C(rowindices);
  transf2C(pivot);
  transf2C(invpivot);

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


