// Backwards Automatic Differentiation of the Cholesky Algorithm
// for calculating partial derivatives of the log-determinant of
// positive definite symmetric sparse matrices.
//
// Implementation by Martin Boer, 2026.
//
// The algorithm combines the sparse Cholesky factorization of
// Ng and Peyton (1993) with the reverse differentiation approach
// of Smith (1995).
//
// References:
//
// Ng, E. G. and B. W. Peyton (1993).
// "Block sparse Cholesky algorithms on advanced uniprocessor computers."
// SIAM Journal on Scientific Computing 14, 1034-1056.
//
// Furrer, R. and S. R. Sain (2010).
// "spam: A sparse matrix R package with emphasis on MCMC
// methods for Gaussian Markov random fields."
// Journal of Statistical Software 36, 1-25.
//
// Smith, S. P. (1995).
// "Differentiation of the Cholesky algorithm."
// Journal of Computational and Graphical Statistics 4, 134-147.
//
// Smith, S. P. (2000).
// "A Tutorial on Simplicity and Computational Differentiation
// for Statisticians."

#include <Rcpp.h>
#include <set>
#include <vector>
#include "AuxFun.h"
#include "SparseMatrix.h"
#include "cholesky.h"
#include "ADcholesky.h"

using namespace Rcpp;
using namespace std;

// helper function for diagXCinvXt
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

