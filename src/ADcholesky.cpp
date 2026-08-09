// Backwards Automated Differentiation of Cholesky Algorithm
// to calculate the partial derivatives of log-determinant of
// positive definite symmetric sparse matrices.
//
// References:
// Ng, Esmond G., and Barry W. Peyton.,
// "Block sparse Cholesky algorithms on advanced uniprocessor computers."
// SIAM Journal on Scientific Computing 14, no. 5 (1993): 1034-1056.
//
// Furrer, Reinhard, and Stephan R. Sain.
// "spam: A sparse matrix R package with emphasis on MCMC
// methods for Gaussian Markov random fields."
// Journal of Statistical Software 36 (2010): 1-25.
//
// Smith, Stephen P. "Differentiation of the Cholesky algorithm."
// Journal of Computational and Graphical Statistics 4, no. 2 (1995): 134-147.
//
// S.P. Smith 2000, A TUTORIAL ON SIMPLICITY AND COMPUTATIONAL DIFFERENTIATION FOR
// STATISTICIANS
//

#include <Rcpp.h>
#include <set>
#include <vector>
#include "AuxFun.h"
#include "SparseMatrix.h"
#include "cholesky.h"

using namespace Rcpp;
using namespace std;

// j is current column in Supernode J
void ADcmod1(NumericVector& F,
             const NumericVector& L, int j, int J,
             const IntegerVector& supernodes,
             const IntegerVector& colpointers)
{
  int s = colpointers[j];
  int e = colpointers[j+1];
  // for all columns in supernode J left to j:
  for (int k=supernodes[J];k<j;k++)
  {
    int jk = colpointers[k] + (j-k);
    int ik = jk;
    double& Fjk = F[jk];
    const double& Ljk = L[jk];
    for (int ij=s; ij<e; ij++)
    {
      // F[ik] = F[ik] - F[ij]*L[jk];
      // F[jk] = F[jk] - F[ij]*L[ik];
      F[ik] -= F[ij]*Ljk;
      Fjk   -= F[ij]*L[ik];
      ik++;
    }
  }
}

// Adjust column j for all columns in supernode K:
void ADcmod2(NumericVector& F,
            const NumericVector& L, int j, int K, int sz,
           NumericVector& t,
           const IntegerVector& indmap,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  // t is dense version of L[j], updated values at end of function:
  int i=0;
  for (int r = rowpointers[K+1] - 1;r>=rowpointers[K];r--)
  {
    int ndx = rowindices[r];
    int pos = colpointers[j+1] - 1 - indmap[ndx];
    t[i++] = F[pos];
  }

  // for all columns k in supernode K:
  for (int k=supernodes[K]; k<supernodes[K+1]; k++)
  {
    int jk = colpointers[k+1]-sz;
    int ik = jk;
    const double& Ljk = L[jk];
    double& Fjk = F[jk];
    for (int i=sz-1;i>=0;i--)
    {
      // F[ik] = F[ik] - F_ij*L[jk];
      // F[jk] = F[jk] - F_ij*L[ik];
      double F_ij = t[i];
      F[ik] -= F_ij*Ljk;
      Fjk   -= F_ij*L[ik];
      ik++;
    }
  }
}

void ADcdiv(NumericVector& F,
            const NumericVector& L, int j, const IntegerVector& colpointers)
{
  const int s = colpointers[j];
  const int e = colpointers[j+1];

  // update AD for column j:
  const double& Ls = L[s];
  double& Fs = F[s];
  for (int i = s + 1; i < e; i++)
  {
    // F[i] = F[i]/L[s];
    // F[s] = F[s] - L[i]*F[i];
    F[i] /= Ls;
    Fs -= L[i]*F[i];
  }
  //F[s] = Fs;
  F[s] = 0.5*F[s]/Ls;
}

void ADcholesky(NumericVector& F,
              const NumericVector& L,
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
  IntegerVector coltop = clone(rowpointers);
  for (int J=0; J<Nsupernodes;J++)
  {
    int szNode = supernodes[J+1] - supernodes[J];
    coltop[J] += szNode-1;
    colhead[J] = rowpointers[J+1]-1;
    if (colhead[J] > coltop[J])
    {
      int rNdx = rowindices[colhead[J]];
      insert(HEAD, LINK, rNdx, J);
    }
  }
  IntegerVector indmap(N,0);
  NumericVector t(N);
  for (int J=Nsupernodes-1; J>=0;J--)
  {
    makeIndMap(indmap, J, rowpointers, rowindices);
    for (int j = supernodes[J+1]-1; j>=supernodes[J]; j--)
    {
      ADcdiv(F, L, j, colpointers);
      ADcmod1(F, L, j, J, supernodes, colpointers);

      int K = HEAD[j];
      while (K!=-1)
      {
        int nextK = LINK[K];
        colhead[K]--;
        if (colhead[K] > coltop[K])
        {
           int rNdx = rowindices[colhead[K]];
           insert(HEAD, LINK, rNdx, K);
        }
        int sz = rowpointers[K+1] - 1 - colhead[K];
        ADcmod2(F, L, j, K, sz, t, indmap, supernodes, rowpointers,colpointers,rowindices);
        K = nextK;
      }
      HEAD[j] = -1;
    }
  }
  return;
}

void initAD(NumericVector& F, const NumericVector& L, const IntegerVector& colpointers)
{
  std::fill(F.begin(), F.end(), 0.0);
  const int N = colpointers.size() - 1;
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    F[s] = 2.0/L[s];
  }
}

void updateH(NumericVector& H, const SparseMatrix& tX, int i, int j, double alpha)
{
  int s1 = tX.rowpointers[i];
  int e1 = tX.rowpointers[i+1];

  int s2 = tX.rowpointers[j];
  int e2 = tX.rowpointers[j+1];

  while (s1 != e1 && s2 != e2)
  {
    if (tX.colindices[s1] < tX.colindices[s2]) ++s1;
    else if (tX.colindices[s1] > tX.colindices[s2]) ++s2;
    else {
      int ndx = tX.colindices[s1]; // = tD.colindices[2]
      H[ndx] += tX.entries[s1]*tX.entries[s2]*alpha;
      ++s1;
      ++s2;
    }
  }
}

// [[Rcpp::export]]
NumericVector diagXCinvXt(Rcpp::S4 obj, Rcpp::S4 transposeX)
{
  SparseMatrix tX(transposeX);
  const int nPred = tX.dim[1];

  IntegerVector supernodes = GetIntVector(obj, "supernodes", 0);
  // Exchange row and columns compared to spam object, as in Ng and Peyton 1993
  IntegerVector colpointers = GetIntVector(obj, "rowpointers", 0);
  IntegerVector rowpointers = GetIntVector(obj, "colpointers", 0);
  IntegerVector rowindices = GetIntVector(obj, "colindices", 0);

  NumericVector L = Rcpp::clone<Rcpp::NumericVector>(obj.slot("entries"));

  const int sz = L.size();
  NumericVector F(sz, 0.0);
  initAD(F, L, colpointers);
  ADcholesky(F, L, supernodes, rowpointers, colpointers, rowindices);

  NumericVector H(nPred, 0.0);

  const int Nsupernodes = supernodes.size()-1;
  for (int J=0; J<Nsupernodes;J++)
  {
    int s = rowpointers[J];
    for (int j=supernodes[J]; j<supernodes[J+1]; j++)
    {
      int k = s;
      for (int ndx = colpointers[j]; ndx < colpointers[j+1]; ndx++)
      {
        int i = rowindices[k++];
        double alpha = F[ndx];
        updateH(H, tX, i, j, alpha);
      }
      s++;
    }
  }
  return H;
}

// [[Rcpp::export]]
List update_Rcpp_fun(Rcpp::S4 obj) {
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector rowindices = obj.slot("rowindices");

  NumericVector L = obj.slot("entries");
  NumericVector F = obj.slot("ADentries");

  cholesky(L, supernodes, rowpointers, colpointers, rowindices);
  initAD(F, L, colpointers);
  ADcholesky(F, L, supernodes, rowpointers, colpointers, rowindices);

  List L_obj;
  L_obj["entries"] = L;
  L_obj["ADentries"] = F;
  return L_obj;
}

// [[Rcpp::export]]
NumericVector solve_Rcpp_fun(Rcpp::S4 obj, const NumericVector& b) {
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericVector L = obj.slot("entries");
  IntegerVector pivot = obj.slot("pivot");
  IntegerVector invpivot = obj.slot("invpivot");

  NumericVector z= forwardCholesky(L, b, supernodes, rowpointers,
                                     colpointers, rowindices, pivot, invpivot);
  NumericVector x = backwardCholesky(L, z, supernodes, rowpointers,
                                       colpointers, rowindices, pivot, invpivot);
  return x;
}

// [[Rcpp::export]]
double logdet_Rcpp_fun(Rcpp::S4 obj) {
  IntegerVector colpointers = obj.slot("colpointers");
  NumericVector L = obj.slot("entries");
  return logdet(L, colpointers);
}


// [[Rcpp::export]]
List constructor_LMMsolver_chol(Rcpp::S4 obj_spam) {
  IntegerVector supernodes = GetIntVector(obj_spam, "supernodes", 0);

  // Exchange row and columns compared to spam object, as in Ng and Peyton 1993
  IntegerVector colpointers = GetIntVector(obj_spam, "rowpointers", 0);
  IntegerVector rowpointers = GetIntVector(obj_spam, "colpointers", 0);
  IntegerVector rowindices = GetIntVector(obj_spam, "colindices", 0);

  IntegerVector pivot = GetIntVector(obj_spam, "pivot", 0);
  IntegerVector invpivot = GetIntVector(obj_spam, "invpivot", 0);
  IntegerVector Dim = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("dimension"));

  NumericVector entries = Rcpp::clone<Rcpp::NumericVector>(obj_spam.slot("entries"));
  const int N_entries = entries.size();

  NumericVector ADentries(N_entries);

  const NumericVector& L = entries;
  NumericVector& F = ADentries;

  initAD(F, L, colpointers);
  ADcholesky(F, L, supernodes, rowpointers, colpointers, rowindices);

  List L_obj;
  L_obj["supernodes"] = supernodes;
  L_obj["colpointers"] = colpointers;
  L_obj["rowpointers"] = rowpointers;
  L_obj["rowindices"] =  rowindices;
  L_obj["pivot"] = pivot;
  L_obj["invpivot"] = invpivot;
  L_obj["entries"] = entries;
  L_obj["ADentries"] = ADentries;
  return L_obj;
}

// Convert a SparseMatrix to the internal AD-Cholesky ordering.
NumericVector convertSparseMatrix(const SparseMatrix& A,
                                  const IntegerVector& supernodes,
                                  const IntegerVector& rowpointers,
                                  const IntegerVector& colpointers,
                                  const IntegerVector& rowindices)
{
  const int Nsupernodes = supernodes.size() - 1;
  const int N = colpointers.size() - 1;
  const int size = colpointers[N];

  NumericVector result(size, 0.0);

  for (int J = 0; J < Nsupernodes; J++)
  {
    for (int j = supernodes[J]; j < supernodes[J + 1]; j++)
    {
      int k   = rowpointers[J + 1] - 1;
      int ndx = colpointers[j + 1] - 1;

      for (int ll = A.rowpointers[j + 1] - 1; ll >= A.rowpointers[j]; ll--)
      {
        int c = A.colindices[ll];

        if (c < j)
          break;

        while (rowindices[k] != c)
        {
          k--;
          ndx--;
        }

        if (k < 0)
        {
          Rcpp::Rcout << "\nPattern mismatch\n";
          Rcpp::Rcout << "Column j = " << j
                      << ", searching for row c = " << c << "\n";

          Rcpp::Rcout << "\nSparseMatrix column rows: ";
          for (int t = A.rowpointers[j]; t < A.rowpointers[j + 1]; t++)
            Rcpp::Rcout << A.colindices[t] << " ";

          Rcpp::Rcout << "\nAD column rows: ";
          for (int t = colpointers[j]; t < colpointers[j + 1]; t++)
            Rcpp::Rcout << rowindices[t] << " ";

          Rcpp::Rcout << "\n";

          Rcpp::stop("Pattern mismatch");
        }

        result[ndx] = A.entries[ll];

        if (c == j)
          break;
      }
    }
  }

  return result;
}

// [[Rcpp::export]]
NumericVector vec(Rcpp::S4 ADobj,
                  Rcpp::S4 spam_matrix)
{
  IntegerVector supernodes = ADobj.slot("supernodes");
  IntegerVector rowpointers = ADobj.slot("rowpointers");
  IntegerVector colpointers = ADobj.slot("colpointers");
  IntegerVector rowindices  = ADobj.slot("rowindices");
  IntegerVector pivot       = ADobj.slot("pivot");

  SparseMatrix A(spam_matrix);

  SparseMatrix Aperm = permuteSymmetric(A,pivot);

  return convertSparseMatrix(
    Aperm,
    supernodes,
    rowpointers,
    colpointers,
    rowindices);
}

