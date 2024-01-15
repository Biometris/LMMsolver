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

using namespace Rcpp;
using namespace std;

// U is a cholesky matrix
// ZtZ is crossproduct design matrix Z
// P is a precision matrix.
// [[Rcpp::export]]
List construct_ADchol_Rcpp(SEXP U,
                           const List& P_list) {
  Rcpp::S4 obj_spam(U);
  IntegerVector supernodes = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("supernodes"));

  // Exchange row and columns compared to spam object, as in Ng and Peyton 1993
  IntegerVector colpointers = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("rowpointers"));
  IntegerVector rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("colpointers"));
  IntegerVector rowindices = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("colindices"));

  IntegerVector pivot = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("pivot"));
  IntegerVector invpivot = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("invpivot"));

  IntegerVector Dim = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("dimension"));

  // copy not really needed or used ...
  NumericVector entries = Rcpp::clone<Rcpp::NumericVector>(obj_spam.slot("entries"));
  NumericVector ADentries = Rcpp::clone<Rcpp::NumericVector>(obj_spam.slot("entries"));

  transf2C(supernodes);
  transf2C(rowpointers);
  transf2C(colpointers);
  transf2C(rowindices);
  transf2C(pivot);
  transf2C(invpivot);

  const int Nsupernodes = supernodes.size()-1;
  const int N = colpointers.size() - 1;
  const int size = colpointers[N];
  const int n_prec_matrices = P_list.size();
  NumericMatrix P_matrix(size, n_prec_matrices);
  for (int i=0;i<n_prec_matrices;i++)
  {
    Rcpp::S4 obj(P_list[i]);
    IntegerVector rowpointers_P = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("rowpointers"));
    IntegerVector colindices_P  = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colindices"));
    NumericVector entries_P     = obj.slot("entries");

    transf2C(rowpointers_P);
    transf2C(colindices_P);

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
  // get number of elements r >= j in supernode K:

  // t is dense version of L[j], updated values at end of function:
  //const int N = colpointers.size() - 1;
  int i=0;
  for (int r = rowpointers[K+1] - 1;r>=rowpointers[K];r--)
  {
    int ndx = rowindices[r];
    int pos = colpointers[j+1] - 1 - indmap[ndx];
    t[i++] = F[pos];
    //i++;
    //if (ndx==j)
    //{
    //  break;
    //}
  }
  //Rcout << "test " << i << " " << sz << endl;
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


void ADcholesky(NumericVector& F,
              const NumericVector& L,
              const IntegerVector& supernodes,
              const IntegerVector& rowpointers,
              const IntegerVector& colpointers,
              const IntegerVector& rowindices)
{
  const int N = colpointers.size() - 1;
  const int Nsupernodes = supernodes.size()-1;
  vector<Node*> S(N);

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
      Node *nodeJ = new Node(J);
      S[rNdx] = add(nodeJ, S[rNdx]);
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
       for (Node **ptr = &S[j]; *ptr;)
       {
          Node *f = removefirstnode(ptr);
          int K = f->data;
          colhead[K]--;
          if (colhead[K] > coltop[K])
          {
             int rNdx = rowindices[colhead[K]];
             S[rNdx] = add(f, S[rNdx]);
          }  else {
             delete f;
          }
          int sz = rowpointers[K+1] - 1 - colhead[K];
          ADcmod2(F, L, j, K, sz, t, indmap, supernodes, rowpointers,colpointers,rowindices);
       }
    }
  }
  return;
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

void initAD(NumericVector& F, const NumericVector& L, const IntegerVector& colpointers)
{
  const int N = colpointers.size() - 1;
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    F[s] = 2.0/L[s];
  }
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


/*
// [[Rcpp::export]]
double logdet(SEXP arg, NumericVector lambda)
{
  Rcpp::S4 obj(arg);
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
  //Rcpp::S4 obj(ADobj);
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericVector L = obj.slot("entries");
  NumericVector F = obj.slot("ADentries");
  NumericMatrix P = obj.slot("P");
  //NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

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
NumericVector partialDerivCholesky(SEXP cholC)
{
  Rcpp::S4 obj(cholC);

  // We use transpose for calculating Automated Differentiation.
  IntegerVector supernodes = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("supernodes"));
  IntegerVector colpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("rowpointers"));
  IntegerVector rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colpointers"));
  IntegerVector rowindices = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colindices"));
  NumericVector L = Rcpp::clone<Rcpp::NumericVector>(obj.slot("entries"));

  // C using indices starting at 0:
  transf2C(supernodes);
  transf2C(colpointers);
  transf2C(rowpointers);
  transf2C(rowindices);

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
NumericVector diagXCinvXt(SEXP cholC, SEXP transposeX)
{
  Rcpp::S4 obj(cholC);
  SparseMatrix tX(transposeX);
  const int nPred = tX.dim[1];

  // We use transpose for calculating Automated Differentiation.
  IntegerVector supernodes = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("supernodes"));
  IntegerVector colpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("rowpointers"));
  IntegerVector rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colpointers"));
  IntegerVector rowindices = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colindices"));
  NumericVector L = Rcpp::clone<Rcpp::NumericVector>(obj.slot("entries"));

  // C using indices starting at 0:
  transf2C(supernodes);
  transf2C(colpointers);
  transf2C(rowpointers);
  transf2C(rowindices);

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


