// Martin Boer, 09 jan 2022
// Biometris, Wageningen University and Research The Netherlands
//
// Backwards Automated Differentation of Cholesky Algorithm
// to calculate the partial derivatives of log-determinant of
// a positive definite symmetric (sparse) matrix.
//
// For details on the implementation of sparse Cholesky, see:
// Ng and Peyton 1993, Furrer and Sain 2010
//
// For details on Backwards Automated Differentiation see:
// S.P. Smith 1995, Differentiation of the Cholesky Algorithm
// S.P. Smith 2000, A TUTORIAL ON SIMPLICITY AND COMPUTATIONAL DIFFERENTIATION FOR
// STATISTICIANS
//
#include <Rcpp.h>
#include <set>
#include <vector>
#include "DiffCholeksyAuxFun.h"

using namespace Rcpp;
using namespace std;

// U is a cholesky matrix
// ZtZ is crossproduct design matrix Z
// P is a precision matrix.
// [[Rcpp::export]]
List construct_ADchol_Rcpp_NgPeyton(SEXP U,
                           const List& P_list) {
  Rcpp::S4 obj_spam(U);
  IntegerVector supernodes = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("supernodes"));
  IntegerVector rowpointers = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("rowpointers"));
  IntegerVector colpointers = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("colpointers"));
  IntegerVector colindices = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("colindices"));
  IntegerVector pivot = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("pivot"));
  IntegerVector Dim = Rcpp::clone<Rcpp::IntegerVector>(obj_spam.slot("dimension"));
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
  //int nelem = Dim[1];
  //Rcout << "Number of precision matrices: " << n_prec_matrices << endl;
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
  L["P"] = P_matrix;
  return L;
}

// make indmap for supernode J:
void makeIndMap(IntegerVector& indmap,
                       int J,
                       const IntegerVector& rowpointers,
                       const IntegerVector& rowindices)
{
  // make indmap for current supernode J:
  int s = rowpointers[J];
  int e  = rowpointers[J+1];
  int l = 0;
  for (int i=e-1; i>=s;i--)
  {
    indmap[rowindices[i]] = l++;
  }
}

// j is current column in Supernode J
void cmod1(NumericVector& L, int j, int J, const IntegerVector& supernodes,
                   const IntegerVector& colpointers)
{
  //if (sj >= J) return 0;
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

void cdiv(NumericVector& L, int j, const IntegerVector& colpointers)
{
  const int s = colpointers[j];
  const int e = colpointers[j+1];
  // pivot:
  L[s] = sqrt(L[s]);
  // update column j:
  for (int i = s + 1; i < e; i++)
  {
    L[i] /= L[s];
  }
}

// clear Set S for columns in supernode J:
void clearSet(vector< set<int> >& S, int J, const IntegerVector& supernodes)
{
  for (int j=supernodes[J];j<supernodes[J+1];j++)
  {
    S[j].clear();
  }
}

vector< set<int> > cholesky(NumericVector& L,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  const int N = colpointers.size() - 1;
  const int Nsupernodes = supernodes.size()-1;

  IntegerVector colhead = clone(rowpointers);
  for (int J=0; J<Nsupernodes;J++)
  {
    int szNode = supernodes[J+1] - supernodes[J];
    colhead[J] += szNode-1;
  }

  IntegerVector indmap(N,0);
  NumericVector t(N);
  vector< set<int> > S(N);
  for (int J=0; J<Nsupernodes;J++) {
    makeIndMap(indmap, J, rowpointers, rowindices);
    for (int j=supernodes[J];j<supernodes[J+1];j++)
    {
      for (set<int>::const_iterator it=S[j].begin(); it!=S[j].end(); it++)
      {
        int K = *it;
        cmod2(L, j, K, t, indmap, supernodes, rowpointers, colpointers, rowindices);

        colhead[K]++;
        if (colhead[K] < rowpointers[K+1]) {
          int rNdx = rowindices[colhead[K]];
          S[rNdx].insert(K);
        }
      }
      cmod1(L, j, J, supernodes, colpointers);
      cdiv(L, j, colpointers);
    }
    colhead[J]++;
    if (colhead[J] < rowpointers[J+1]) {
      int rNdx = rowindices[colhead[J]];
      S[rNdx].insert(J);
    }
    //clearSet(S, J, supernodes);
  }
  return S;
}


void ADcholesky(NumericVector& F,
              const NumericVector& L,
              vector<set<int> >& S,
              const IntegerVector& supernodes,
              const IntegerVector& rowpointers,
              const IntegerVector& colpointers,
              const IntegerVector& rowindices)
{
  cout << "outline ADcholeksy" << endl;
  const int N = colpointers.size() - 1;
  const int Nsupernodes = supernodes.size()-1;

  IntegerVector colhead = clone(rowpointers);
  IntegerVector coltop = clone(rowpointers);
  for (int J=0; J<Nsupernodes;J++)
  {
    int szNode = supernodes[J+1] - supernodes[J];
    coltop[J] += szNode-1;
    colhead[J] = rowpointers[J+1] - 1;
  }
  IntegerVector indmap(N,0);
  NumericVector t(N);
  for (int J=Nsupernodes-1; J>=0;J--)
  {
    cout << " supernode " << J << endl;
    makeIndMap(indmap, J, rowpointers, rowindices);
    for (int j = supernodes[J+1]-1; j>=supernodes[J]; j--)
    {
       cout << "  column " << j << endl;
       cout << "    ADcdiv; ADcmod1" << endl;
       for (set<int>::const_iterator it=S[j].begin(); it!=S[j].end(); it++)
       {
         cout << "    ADcmod2 supernode " << *it << endl;
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

NumericVector initAD(const NumericVector& L, const IntegerVector& colpointers)
{
  const int N = colpointers.size() - 1;
  const int sz = L.size();
  NumericVector F(sz, 0.0);
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    F[s] = 2.0/L[s];
  }
  return F;
}


// [[Rcpp::export]]
double logdetNgPeyton(SEXP arg, NumericVector lambda)
{
  Rcpp::S4 obj(arg);
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

  // define matrix L (lower triangle matrix values)
  const int sz = P.nrow();
  const int n_prec_mat = P.ncol();
  NumericVector L(sz, 0.0);
  for (int k=0;k<n_prec_mat;k++)
  {
    NumericVector Pk = P(_, k);
    double alpha = lambda[k];
    for (int i=0;i<sz;i++)
    {
      L[i] += alpha*Pk[i];
    }
  }
  vector<set<int> > S = cholesky(L, supernodes, rowpointers, colpointers, rowindices);
  return logdet(L, colpointers);
}


// [[Rcpp::export]]
NumericVector dlogdetNgPeyton(SEXP arg, NumericVector lambda)
{
  Rcpp::S4 obj(arg);
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

  // define matrix L (lower triangle matrix values)
  const int sz = P.nrow();
  const int n_prec_mat = P.ncol();
  NumericVector L(sz, 0.0);
  for (int k=0;k<n_prec_mat;k++)
  {
    NumericVector Pk = P(_, k);
    double alpha = lambda[k];
    for (int i=0;i<sz;i++)
    {
      L[i] += alpha*Pk[i];
    }
  }
  vector<set<int> > S = cholesky(L, supernodes, rowpointers, colpointers, rowindices);
  double logDet = logdet(L, colpointers);

  NumericVector F = initAD(L, colpointers);
  ADcholesky(F, L, S, supernodes, rowpointers, colpointers, rowindices);


  // calculate the partial derivatives:
  const int N = colpointers.size()-1;
  NumericVector gradient(n_prec_mat);
  for (int i=0;i<F.size();i++)
  {
    for (int k=0;k<n_prec_mat;k++)
      gradient[k] += F[i]*P(i,k);
  }
  double sum = 0.0;
  for (int i=0;i<n_prec_mat;i++)
  {
    sum += lambda[i]*gradient[i];
  }
  // correction, to make inproduct
  // between lambda and gradient equal to N:
  for (int i=0;i<n_prec_mat;i++)
  {
    gradient[i] *= N/sum;
  }

  //NumericVector gradient(n_prec_mat,0.0);
  gradient.attr("logdet") = logDet;

  return gradient;
}





