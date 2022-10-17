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
  IntegerVector rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("rowpointers"));
  IntegerVector colpointers = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("colpointers"));
  IntegerVector colindices = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("colindices"));
  IntegerVector pivot = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("pivot"));
  IntegerVector Dim = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("dimension"));

  // copy not really needed or used ...
  NumericVector entries = Rcpp::clone<Rcpp::NumericVector>(obj_spam.slot("entries"));
  NumericVector ADentries = Rcpp::clone<Rcpp::NumericVector>(obj_spam.slot("entries"));

  transf2C(supernodes);
  transf2C(rowpointers);
  transf2C(colpointers);
  transf2C(colindices);
  transf2C(pivot);

  const int size = rowpointers[rowpointers.size()-1];

  const int Nnodes = supernodes.size() - 1;
  std::vector<int> colindices_ext(size);
  int cnt = 0;
  for (int i=0;i<Nnodes;i++)
  {
    int sz_node = supernodes[i+1] - supernodes[i];
    int s_col = colpointers[i];
    int e_col = colpointers[i+1];
    for (int k=0;k<sz_node;k++)
    {
      for (int l=s_col+k;l<e_col;l++)
      {
        colindices_ext[cnt++] = colindices[l];
      }
    }
  }

  cnt = 0;
  std::vector<int> rowindices_ext(size);
  const int dim = rowpointers.size() - 1;
  for (int i=0;i<dim;i++)
  {
    int n_elem_cur = rowpointers[i+1] - rowpointers[i];
    for (int k=0;k<n_elem_cur;k++)
    {
      rowindices_ext[cnt++] = i;
    }
  }

  const int n_prec_matrices = P_list.size();
  NumericMatrix P_matrix(size, n_prec_matrices);
  for (int i=0;i<n_prec_matrices;i++)
  {
    Rcpp::S4 obj_P(P_list[i]);
    vector<double> entries_P = convert_matrix(rowindices_ext, colindices_ext, pivot, obj_P);
    for (int j=0;j<size;j++)
    {
      P_matrix(j,i) = entries_P[j];
    }
  }

  List L;
  L["supernodes"] = supernodes;
  L["colpointers"] = rowpointers;
  L["rowpointers"] = colpointers;
  L["rowindices"] =  colindices;
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
  int e  = rowpointers[J+1];
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
  int s = colpointers[j];
  int e = colpointers[j+1];
  // for all columns in supernode J left to j:
  for (int k=supernodes[J];k<j;k++)
  {
    int jk = colpointers[k] + (j-k);
    int ik = jk;
    double Ljk = L[jk];
    for (int ij=s; ij<e; ij++)
    {
       L[ij] -= L[ik]*Ljk;
       ik++;
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
void cmod2(NumericVector& L, int j, int K,
           NumericVector& t,
           const IntegerVector& indmap,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  // get number of elements r >= j in supernode K:
  int sz =0;

  // t is dense version of L[j], updated values at end of function:
  //const int N = colpointers.size() - 1;
  for (int r = rowpointers[K+1] - 1;r>=rowpointers[K];r--)
  {
    int ndx = rowindices[r];
    //int pos = colpointers[j+1] - 1 - indmap[ndx];
    t[sz] = 0.0;
    sz++;
    if (ndx==j)
    {
      break;
    }
  }

  // for all columns k in supernode K:
  for (int k=supernodes[K]; k<supernodes[K+1]; k++)
  {
    int jk = colpointers[k+1]-sz;
    int ik = jk;
    double Ljk = L[jk];
    for (int i=sz-1;i>=0;i--)
    {
      t[i] += L[ik]*Ljk;
      ik++;
    }
  }

  // write results dense matrix t back to L_j
  sz=0;
  for (int r = rowpointers[K+1] - 1;r>=rowpointers[K];r--)
  {
    int ndx = rowindices[r];
    int pos = colpointers[j+1] - 1 - indmap[ndx];
    L[pos] -= t[sz];
    sz++;
    if (ndx==j)
    {
      break;
    }
  }
}

// Adjust column j for all columns in supernode K:
void ADcmod2(NumericVector& F,
            const NumericVector& L, int j, int K,
           NumericVector& t,
           const IntegerVector& indmap,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  // get number of elements r >= j in supernode K:
  int sz =0;

  // t is dense version of L[j], updated values at end of function:
  //const int N = colpointers.size() - 1;
  for (int r = rowpointers[K+1] - 1;r>=rowpointers[K];r--)
  {
    int ndx = rowindices[r];
    int pos = colpointers[j+1] - 1 - indmap[ndx];
    t[sz] = F[pos];
    sz++;
    if (ndx==j)
    {
      break;
    }
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

void cdiv(NumericVector& L, int j, const IntegerVector& colpointers)
{
  const int s = colpointers[j];
  const int e = colpointers[j+1];
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
        cmod2(L, j, K, t, indmap, supernodes, rowpointers, colpointers, rowindices);

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
          ADcmod2(F, L, j, K, t, indmap, supernodes, rowpointers,colpointers,rowindices);
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
NumericVector dlogdet(SEXP ADobj, NumericVector theta)
{
  Rcpp::S4 obj(ADobj);
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

/*

Not used, just to show the structure of the sparse cholesky matrix with supernodes.

// [[Rcpp::export]]
NumericMatrix PrintCholesky(SEXP cholC)
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
