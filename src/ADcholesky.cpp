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

// U is a cholesky matrix
// ZtZ is crossproduct design matrix Z
// P is a precision matrix.
// [[Rcpp::export]]
List construct_ADchol_Rcpp(Rcpp::S4 obj_spam,
                           const List& P_list) {
  //Rcpp::S4 obj_spam(U);
  IntegerVector supernodes = GetIntVector(obj_spam, "supernodes", 0);

  // Exchange row and columns compared to spam object, as in Ng and Peyton 1993
  IntegerVector colpointers = GetIntVector(obj_spam, "rowpointers", 0);
  IntegerVector rowpointers = GetIntVector(obj_spam, "colpointers", 0);
  IntegerVector rowindices = GetIntVector(obj_spam, "colindices", 0);

  IntegerVector pivot = GetIntVector(obj_spam, "pivot", 0);
  IntegerVector invpivot = GetIntVector(obj_spam, "invpivot", 0);

  IntegerVector Dim = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("dimension"));

  // copy not really needed or used ...
  NumericVector entries = Rcpp::clone<Rcpp::NumericVector>(obj_spam.slot("entries"));
  NumericVector ADentries = Rcpp::clone<Rcpp::NumericVector>(obj_spam.slot("entries"));

  const int Nsupernodes = supernodes.size()-1;
  const int N = colpointers.size() - 1;
  const int size = colpointers[N];
  const int n_prec_matrices = P_list.size();
  NumericMatrix P_matrix(size, n_prec_matrices);
  for (int i=0;i<n_prec_matrices;i++)
  {
    Rcpp::S4 obj(P_list[i]);
    IntegerVector rowpointers_P = GetIntVector(obj, "rowpointers", 0);
    IntegerVector colindices_P  = GetIntVector(obj, "colindices", 0);
    NumericVector entries_P     = obj.slot("entries");

    vector<double> result(size, 0.0);
    vector<double> z(N);
    for (int J=0; J<Nsupernodes;J++)
    {
      int s = rowpointers[J];
      for (int j=supernodes[J]; j<supernodes[J+1]; j++)
      {
        int r = pivot[j];
        if (rowpointers_P[r] != rowpointers_P[r+1]) {
          std::fill(z.begin(), z.end(), 0.0);
          for (int ll=rowpointers_P[r];ll<rowpointers_P[r+1];ll++) {
            int c = colindices_P[ll];
            z[invpivot[c]] = entries_P[ll];
          }

          int k = s;
          for (int ndx = colpointers[j]; ndx < colpointers[j+1]; ndx++)
          {
            int ii = rowindices[k++];
            result[ndx] = z[ii];
          }
        }
        s++;
      }
    }

    for (int j=0;j<size;j++)
    {
      P_matrix(j,i) = result[j];
    }
  }

  List L;
  L["supernodes"] = supernodes;
  L["colpointers"] = colpointers;
  L["rowpointers"] = rowpointers;
  L["rowindices"] =  rowindices;
  L["pivot"] = pivot;
  L["invpivot"] = invpivot;
  L["entries"] = entries;
  L["ADentries"] = ADentries;
  L["P"] = P_matrix;
  return L;
}

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
  const int N = colpointers.size() - 1;
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    F[s] = 2.0/L[s];
  }
}

//' Calculate the partial derivatives of log-determinant.
//'
//' This function calculates the partial derivatives of the the log-determinant in an
//' efficient way, by using reverse Automated Differentiation of the Cholesky Algorithm,
//' see Smith (1995) for details. Let
//'  \deqn{C = \sum_{i} \theta_i P_i}
//' where the matrices \eqn{P_i} are stored in the `ADchol` object. The partial derivatives
//' of matrix \eqn{C} are defined by:
//' \deqn{\frac{\partial C}{\partial \theta_i} = \text{trace} [C^{-1} P_i]},
//' but are calculated in a more efficient way using backwards Automated Differentiation.
//'
//' @param ADobj object of class ADchol.
//' @param theta a vector with precision or penalty parameters
//'
//' @return The gradient with partial derivatives of \eqn{log|C|} with respect to
//' parameters \eqn{\theta_i}. As attribute \code{logdet}, \eqn{log|C|} is returned.
//'
//' @references
//' Smith, S. P. (1995). Differentiation of the Cholesky algorithm.
//' Journal of Computational and Graphical Statistics, 4(2), 134-147.
//'
//' @noRd
//' @keywords internal
//'
// [[Rcpp::export]]
NumericVector dlogdet(Rcpp::S4 obj, NumericVector theta,
                      Nullable<NumericVector> b_ = R_NilValue)
{
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericVector L = obj.slot("entries");
  NumericVector F = obj.slot("ADentries");
  NumericMatrix P = obj.slot("P");

  // define matrix L (lower triangle matrix values)
  const int sz = P.nrow();

  const int n_prec_mat = P.ncol();

  if (n_prec_mat != theta.size()) {
    stop("wrong length vector theta ");
  }

  std::fill(L.begin(), L.end(), 0.0);
  std::fill(F.begin(), F.end(), 0.0);

  for (int k=0;k<n_prec_mat;k++)
  {
    NumericMatrix::Column Pk = P(_, k);
    double alpha = theta[k];
    for (int i=0;i<sz;i++)
    {
      L[i] += alpha*Pk[i];
    }
  }
  cholesky(L, supernodes, rowpointers, colpointers, rowindices);
  double logDet = logdet(L, colpointers);

  //NumericVector F(sz, 0.0);
  initAD(F, L, colpointers);
  ADcholesky(F, L, supernodes, rowpointers, colpointers, rowindices);

  NumericVector gradient(n_prec_mat);
  for (int k=0;k<n_prec_mat;k++)
  {
    NumericMatrix::Column Pk = P(_, k);
    gradient[k] = std::inner_product(F.begin(), F.end(), Pk.begin(), 0.0);
  }

  // correction, to make sure inner product
  // between theta and gradient equal to N:
  double sum = 0.0;
  const int N = colpointers.size()-1;
  for (int i=0;i<n_prec_mat;i++)
  {
    sum += theta[i]*gradient[i];
  }
  for (int i=0;i<n_prec_mat;i++)
  {
    gradient[i] *= N/sum;
  }

  gradient.attr("logdet") = logDet;

  if (b_.isNotNull()) {
    NumericVector b(b_);        // casting to underlying type NumericVector
    IntegerVector pivot = obj.slot("pivot");
    IntegerVector invpivot = obj.slot("invpivot");

    NumericVector z= forwardCholesky(L, b, supernodes, rowpointers,
                             colpointers, rowindices, pivot, invpivot);
    NumericVector x = backwardCholesky(L, z, supernodes, rowpointers,
                     colpointers, rowindices, pivot, invpivot);
    gradient.attr("x.coef") = x;
  }
  return gradient;
}

// [[Rcpp::export]]
NumericVector partialDerivCholesky(Rcpp::S4 obj)
{
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
  return F;
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

/*

// [[Rcpp::export]]
NumericVector ForwardCholesky(SEXP cholC, NumericVector& b)
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

  return forwardCholesky(L, b, supernodes, rowpointers,
                      colpointers, rowindices, pivot, invpivot);

}


// [[Rcpp::export]]
NumericVector BackwardCholesky(SEXP cholC, NumericVector& b)
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

  return backwardCholesky(L, b, supernodes, rowpointers,
                         colpointers, rowindices, pivot, invpivot);
}

*/


