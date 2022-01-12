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


vector< set<int> > makeSetS(const IntegerVector& supernodes,
                          const IntegerVector& rowpointers,
                         const IntegerVector& colpointers,
                         const IntegerVector& rowindices)
{
  const int Nsupernodes = supernodes.size()-1;
  const int N = colpointers.size()-1;
  vector< set<int> > S(N);
  Rcout << "make sets S:" << endl;

  // for each supernode:
  for (int s=0; s<Nsupernodes;s++) {
    // for column j:
    for (int j=supernodes[s];j<supernodes[s+1];j++)
    {
      int sz = supernodes[s+1] - supernodes[s];
      for (int r = rowpointers[s]+sz; r < rowpointers[s+1];r++)
      {
        int ndx = rowindices[r];
        S[ndx].insert(s);
      }

      // diagnostics....
      Rcout << "   set  S[" << j << "] = { ";
      for (set<int>::const_iterator it=S[j].begin();it!=S[j].end();it++)
      {
        Rcout << *it << " ";
      }
      Rcout << "}" << endl;
    }
  }
  Rcout << endl;
  return S;
}

vector<int> makeIntMap(int s,
                       const IntegerVector& rowpointers,
                       const IntegerVector& rowindices)
{
  // make indmap for current supernode s:
  int start = rowpointers[s];
  int end  = rowpointers[s+1];
  int sz = end-start;
  vector<int> indmap(sz);
  int l = 0;
  for (int i=end-1; i>=start;i--)
  {
    indmap[l++] = rowindices[i];
  }
  /*
  Rcout << "  Indmap supernode " << s << ": ";
  for (vector<int>::const_iterator it=indmap.begin();it!=indmap.end();it++)
  {
    Rcout << " " << *it;
  }
  Rcout << endl;
  */
  return indmap;
}

// print current supernode and size:
void PrintSuperNodeInfo(int s,
                        const IntegerVector& supernodes,
                        const IntegerVector& rowpointers,
                        const IntegerVector& colpointers,
                        const IntegerVector& rowindices)
{
  Rcout << "SuperNode " << s << " of size "
        << supernodes[s+1] - supernodes[s] << ":" << endl;

  int startRowPtr = rowpointers[s];
  // for column j:
  for (int j=supernodes[s];j<supernodes[s+1];j++)
  {
    Rcout << " column " << j << ": ";
    int r = startRowPtr;
    for (int i=colpointers[j];i<colpointers[j+1];i++)
    {
      Rcout << "r" << rowindices[r] << "[k=" << i << "] ";
      r++;
    }
    Rcout << endl;
    // increase startRowPtr:
    startRowPtr++;
  }
}


// just a very simple test...
void mult2(NumericVector& L, double alpha)
{
  const int N = L.size();
  for (int i=0;i<N;i++) {
    L[i] *= alpha;
  }
}

// [[Rcpp::export]]
double PrintADchol(SEXP arg, NumericVector lambda)
{
  Rcpp::S4 obj(arg);
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

  // give some info:
  Rcout << "P.nrow " << P.nrow() << endl;
  Rcout << "P.ncol " << P.ncol() << endl;
  Rcout << "superNodes:  " <<  supernodes << endl;
  Rcout << "rowpointers: " <<  rowpointers << endl;
  Rcout << "colpointers: " <<  colpointers << endl;
  Rcout << "rowindices:  " <<  rowindices << endl;
  Rcout << endl;

  // define matrix L (lower triangle matrix values)
  const int sz = P.nrow();
  const int n_prec_mat = P.ncol();
  NumericVector L(sz, 0.0);
  for (int i=0;i<sz;i++)
  {
    for (int k=0;k<n_prec_mat;k++)
    {
      L[i] += lambda[k]*P(i,k);
    }
  }

  // simple test.....
  NumericVector L2 = Rcpp::clone<Rcpp::NumericVector>(L);
  mult2(L2, 2.0);
  for (int i=0;i<L.size();i++) {
    Rcout << setw(3) << i << setw(10) << L[i] << setw(10) << L2[i] << endl;
  }

  // make set S_j for each column j, see Ng and Peyton:
  vector<set<int> > S = makeSetS(supernodes, rowpointers, colpointers, rowindices);

  const int Nsupernodes = supernodes.size()-1;
  for (int s=0; s<Nsupernodes;s++) {
    Rcout << "for supernode " << s << endl;
    vector<int> indmap = makeIntMap(s, rowpointers, rowindices);
    for (int j=supernodes[s];j<supernodes[s+1];j++)
    {
      Rcout << "  for column " << j << endl;
      for (set<int>::const_iterator it=S[j].begin(); it!=S[j].end(); it++)
      {
        Rcout << "    correct for columns in supernode " << *it << endl;
      }
      Rcout << "    correct for columns in current supernode " << s << endl;
      Rcout << "    pivot; cdiv" << endl;
    }
  }

  return 0.0;
}


