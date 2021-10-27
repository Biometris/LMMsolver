#include "DiffCholeksy.h"
#include <set>
#include <vector>

using namespace Rcpp;
using namespace std;

// Transform to C++ Notation indices
void transf2C(IntegerVector& ndx)
{
  for (int i=0;i<ndx.size();i++)
  {
    ndx[i] -= 1;
  }
}

// not very efficient (but not too bad): Make a Class?
double getvalueC(IntegerVector rowpointers,
                 IntegerVector colindices,
                 NumericVector entries,
                 int i,
                 int j)
{
  int s = rowpointers[i];
  int e = rowpointers[i+1];
  for (int k=s;k<e;k++)
  {
    if (colindices[k] == j) return entries[k];
  }
  return 0.0;
}

// left looking cholesky
//
// [[Rcpp::export]]
NumericVector cholesky(NumericVector L,
                       const IntegerVector& colpointers,
                       const IntegerVector& rowindices)
{
  const int N = colpointers.size()-1;

  // for details of S_j, see Ng and Peyton, Figure 2.2
  vector<set<int> > S(N);
  // current headings of columns:
  IntegerVector colhead = clone(colpointers);

  IntegerVector indmap(N,0);
  // dense map t,
  vector<double> t(N,0.0);

  for (int j=0;j<N;j++)
  {
    int s = colpointers[j];
    int e = colpointers[j+1];
    int sz = e-s;

    // make indmap for current column j
    // see Ng and Peyton 1993, p. 1040-1041 for details:
    int l = 0;
    for (int i=e-1; i>=s;i--)
    {
      indmap[rowindices[i]] = l++;
    }

    // init dense vector, see Ng and Peyton 1993 for details.
    for (int i=0;i<sz;i++)
    {
      t[i] = 0.0;
    }

    // for each k in S_j:
    for (set<int>::const_iterator it=S[j].begin();it!=S[j].end();it++)
    {
      int k = *it;
      int jk = colhead[k];
      double L_jk = L[jk];
      for (int i = colhead[k]; i < colpointers[k+1];i++)
      {
        int ndx = indmap[rowindices[i]];
        t[ndx] += L[i]*L_jk;
      }
      colhead[k]++;
      if (colhead[k] < colpointers[k+1]) {
        int rNdx = rowindices[colhead[k]];
        S[rNdx].insert(k);
      }
    }

    // update column j with dense vector t:
    int cnt = 0;
    for (int i = e-1; i>=s ;i--)
    {
      L[i] -= t[cnt++];
    }

    // pivot:
    L[s] = sqrt(L[s]);
    // update column j:
    for (int i = s + 1; i < e; i++)
    {
      L[i] /= L[s];
    }
    colhead[j]++;
    if (colhead[j] < colpointers[j+1]) {
      int rNdx = rowindices[colhead[j]];
      S[rNdx].insert(j);
    }
    S[j].clear();
  }
  return L;
}

NumericVector AD_cholesky(const NumericVector& L,
                          const IntegerVector& colpointers,
                          const IntegerVector& rowindices)
{
  const int N = colpointers.size()-1;

  vector<set<int> > S(N);

  IntegerVector indmap(N,0);
  // dense map t,
  vector<double> t(N,0.0);

  // reverse automatic differentiation:
  // set the colheads and sets.
  IntegerVector colhead(N);
  for (int j=0;j<N;j++) {
    colhead[j] = colpointers[j+1]-1;
    int ndx = colhead[j];
    int r = rowindices[ndx];
    if (colhead[j]!=colpointers[j])
    {
      S[r].insert(j);
    }
  }

  // initialize F:
  NumericVector F(L.size());
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    F[s] = 2.0/L[s];
  }

  for (int j=N-1; j>=0; j--)
  {
    int s = colpointers[j];
    int e = colpointers[j+1];

    int l = 0;
    for (int i=e-1; i>=s;i--)
    {
      indmap[rowindices[i]] = l++;
    }

    for (int i = s + 1; i < e; i++)
    {
      //L[i] /= L[s];
      F[i] = F[i]/L[s];
      F[s] = F[s] - L[i]*F[i];
    }
    F[s] = 0.5*F[s]/L[s];

    int cnt=0;
    for (int i = e-1; i>=s;i--)
    {
      t[cnt++] = F[i];
    }

    // for each k in S_j:
    for (set<int>::const_iterator it=S[j].begin();it!=S[j].end();it++)
    {
      int k = *it;
      int jk = colhead[k];
      for (int ik = colhead[k]; ik < colpointers[k+1];ik++)
      {
        int ndx = indmap[rowindices[ik]];
        double F_ij = t[ndx];
        F[ik] = F[ik] - F_ij*L[jk];
        F[jk] = F[jk] - F_ij*L[ik];
      }
      colhead[k]--;
      if (colhead[k] > colpointers[k]) {
        int rNdx = rowindices[colhead[k]];
        S[rNdx].insert(k);
      }
    }
    S[j].clear();
  }
  return F;
}

// Calculate log determinant using Left-looking Cholesky:
//
// [[Rcpp::export]]
double logdet(SEXP arg,
              NumericVector lambda)
{
  Rcpp::S4 obj(arg);

  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

  const int sz = rowindices.size();
  const int n_prec_mat = P.ncol();
  NumericVector C(sz,0.0);

  for (int i=0;i<sz;i++)
  {
    for (int k=0;k<n_prec_mat;k++)
    {
      C[i] += lambda[k]*P(i,k);
    }
  }

  NumericVector L = cholesky(C,colpointers, rowindices);

  const int N = colpointers.size()-1;
  double sum = 0;
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    sum += 2.0*log(L[s]);
  }
  return sum;
}

// backwards Automatic Differentiation, using Left-looking Cholesky:
//
// [[Rcpp::export]]
NumericVector dlogdet(SEXP arg,
                      NumericVector lambda)
{
  Rcpp::S4 obj(arg);

  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

  const int sz = rowindices.size();
  const int n_prec_mat = P.ncol();
  NumericVector C(sz, 0.0);
  for (int i=0;i<sz;i++)
  {
    for (int k=0;k<n_prec_mat;k++)
    {
      C[i] += lambda[k]*P(i,k);
    }
  }

  //Rcout << "start Cholesky" << endl;
  NumericVector L = cholesky(C,colpointers, rowindices);
  //Rcout << "start AD Cholesky" << endl;
  NumericVector F = AD_cholesky(L,colpointers, rowindices);
  //Rcout << "end AD Cholesky" << endl;

  // evaluate dL/dlambda:
  NumericVector gradient(n_prec_mat);
  for (int i=0;i<F.size();i++)
  {
    for (int k=0;k<n_prec_mat;k++)
      gradient[k] += F[i]*P(i,k);
  }

  return gradient;
}


// backwards Automatic Differentiation, using Left-looking Cholesky:
// returns both the logdet and the first derivatives.
//
// [[Rcpp::export]]
NumericVector logdetPlusDeriv(SEXP arg,
                      NumericVector lambda)
{
  Rcpp::S4 obj(arg);

  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

  const int sz = rowindices.size();
  const int n_prec_mat = P.ncol();
  NumericVector C(sz, 0.0);
  for (int i=0;i<sz;i++)
  {
    for (int k=0;k<n_prec_mat;k++)
    {
      C[i] += lambda[k]*P(i,k);
    }
  }

  NumericVector L = cholesky(C, colpointers, rowindices);
  NumericVector F = AD_cholesky(L, colpointers, rowindices);

  // evaluate dL/dlambda:
  NumericVector result(n_prec_mat+1);

  // 1. calculate the logdet
  const int N = colpointers.size()-1;
  double sum = 0;
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    sum += 2.0*log(L[s]);
  }
  result[0] = sum;

  // 2. calculate the derivatives
  NumericVector gradient(n_prec_mat);
  for (int i=0;i<F.size();i++)
  {
    for (int k=0;k<n_prec_mat;k++)
      result[k+1] += F[i]*P(i,k);
  }

  return result;
}


/*
 // forward Automatic Differentation, using Left-looking Cholesky:
 // [[Rcpp::export]]
 double dlogdet_AD_forward(double lambda , SEXP arg)
 {
 Rcpp::S4 obj(arg);

 IntegerVector colpointers = obj.slot("colpointers");
 IntegerVector rowindices = obj.slot("rowindices");
 NumericVector entriesA = Rcpp::clone<Rcpp::NumericVector>(obj.slot("ZtZ"));
 NumericVector entriesB = Rcpp::clone<Rcpp::NumericVector>(obj.slot("P"));

 NumericVector L = entriesA + lambda * entriesB;
 NumericVector F = entriesB;

 const int N = colpointers.size()-1;

 // current headings of columns:
 IntegerVector colhead = clone(colpointers);
 for (int j=0;j<N;j++)
 {
 // see Ng and Peyton 1993, p. 1040-1041 for details:
 map<int,int> indmap = make_indmap(colpointers,rowindices,j);

 int s = colpointers[j];
 int e = colpointers[j+1];
 int sz = e-s;

 vector<double> t0(sz,0.0), t1(sz,0.0);

 for (int k=0;k<j;k++)
 {
 if (rowindices[colhead[k]] == j)
 {
 int jk = colhead[k];
 for (int i = colhead[k]; i < colpointers[k+1];i++)
 {
 int ndx = indmap[rowindices[i]];
 t0[ndx] += L[i]*L[jk];
 t1[ndx] += L[i]*F[jk] + F[i]*L[jk]; // deriv
 }
 colhead[k]++;
 }
 }

 // update column j with dense vector t:
 int cnt = 0;
 for (int i = e-1; i>=s ;i--)
 {
 L[i] -= t0[cnt];
 F[i] -= t1[cnt];
 cnt++;
 }

 // pivot:
 L[s] = sqrt(L[s]);
 F[s] = 0.5*F[s]/L[s];

 // update column j:
 for (int i = s + 1; i < e; i++)
 {
 L[i] /= L[s];
 F[i] = (F[i] - L[i]*F[s])/L[s];
 }
 colhead[j]++;
 }

 double sum = 0;
 for (int k=0;k<N;k++)
 {
 int s = colpointers[k];
 sum += 2.0*F[s]/L[s];
 }
 return sum;
 }

 */

vector<double> convert_matrix(const vector<int>& rowindices_ext,
                              const vector<int>& colindices_ext,
                              const IntegerVector& pivot,
                              SEXP A)
{
  Rcpp::S4 obj(A);
  IntegerVector rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("rowpointers"));
  IntegerVector colindices  = Rcpp::clone<Rcpp::IntegerVector>(obj.slot("colindices"));
  NumericVector entries     = obj.slot("entries") ;

  transf2C(rowpointers);
  transf2C(colindices);

  const int Nrow = rowindices_ext.size();
  std::vector<double> z(Nrow);
  for (int i=0;i<Nrow;i++)
  {
    int r = rowindices_ext[i];
    int c = colindices_ext[i];
    double x = getvalueC(rowpointers,colindices,entries,pivot[r],pivot[c]);
    z[i] = x;
  }
  return z;
}

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
  transf2C(supernodes);
  transf2C(rowpointers);
  transf2C(colpointers);
  transf2C(colindices);
  transf2C(pivot);

  //const int dim =
  const int size = rowpointers[rowpointers.size()-1];
  //cout << "Nelem: " << size  << endl;

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
  //Rcout << "Start rowindices" << endl;

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
  int nelem = rowindices_ext.size();
  //Rcout << "Number of precision matrices: " << n_prec_matrices << endl;
  NumericMatrix P_matrix(nelem, n_prec_matrices);
  for (int i=0;i<n_prec_matrices;i++)
  {
    Rcpp::S4 obj_P(P_list[i]);
    vector<double> entries_P = convert_matrix(rowindices_ext, colindices_ext, pivot, obj_P);
    for (int j=0;j<nelem;j++)
    {
      P_matrix(j,i) = entries_P[j];
    }
  }

  List L;
  L["colpointers"] = rowpointers;
  L["rowindices"] = Rcpp::wrap(colindices_ext);
  L["P"] = P_matrix;
  return L;
}
