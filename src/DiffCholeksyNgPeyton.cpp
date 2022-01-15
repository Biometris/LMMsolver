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

// make intmap for supernode J:
map<int, int> makeIntMap(int J, int N,
                       const IntegerVector& rowpointers,
                       const IntegerVector& rowindices)
{
  // make indmap for current supernode s:
  int s = rowpointers[J];
  int e  = rowpointers[J+1];
  //int sz = end-start;
  map<int, int> indmap;
  int l = 0;
  for (int i=e-1; i>=s;i--)
  {
    indmap[rowindices[i]] = l++;
  }

  //int l = 0;
  //for (int i=end-1; i>=start;i--)
  //{
  //  indmap[l++] = rowindices[i];
  //}

  Rcout << "  Indmap supernode " << J << ": ";
  for (map<int,int>::const_iterator it=indmap.begin();it!=indmap.end();it++)
  {
    Rcout << " [k=" << it->first << ",v=" << it->second << "] ";
  }
  Rcout << endl;

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

// j is current column in Supernode J
void Cmod_inside_J(NumericVector& L, int j, int J,
                   const IntegerVector& supernodes,
                   const IntegerVector& colpointers)
{
  Rcout << "    Cmod_inside_J: correct for columns in current supernode " << J << endl;

  // first simple check, without calculations:
  // return;

  // for all columns in supernode J left to j:
  for (int k=supernodes[J];k<j;k++)
  {
    int jk = colpointers[k] + (j-k);
    int ik = jk;
    for (int ij=colpointers[j]; ij<colpointers[j+1]; ij++)
    {
       L[ij] = L[ij] - L[ik]*L[jk];
       ik++;
    }
  }
}

// Adjust column j for all columns in supernode K:
void Cmod_outside_J(NumericVector& L, int j, int K,
                    map<int,int>& indmap,
                            const IntegerVector& supernodes,
                            const IntegerVector& rowpointers,
                            const IntegerVector& colpointers,
                            const IntegerVector& rowindices)
{
  Rcout << "    Cmod_outside_J: correct for columns in supernode " << K << endl;
  Rcout << "      Indmap supernode: ";
  for (map<int,int>::const_iterator it=indmap.begin();it!=indmap.end();it++)
  {
    Rcout << " [k=" << it->first << ",v=" << it->second << "] ";
  }
  Rcout << endl;

  // get number of elements r >= j in supernode K:
  int sz =0;

  // t is dense version of L[j], updated values at end of function:
  NumericVector t;
  IntegerVector posL;
  //Rcout << "      init dense vector t " << endl;
  for (int r = rowpointers[K+1] - 1;r>=rowpointers[K];r--)
  {
    int ndx = rowindices[r];

    //int e = colpointers[j+1];
    int pos = colpointers[j+1] - 1 - indmap[ndx];
    posL.push_back(pos);
    t.push_back(L[pos]);

    Rcout << "         in init dense:" << setw(3) << r
          << setw(3) << ndx << setw(3) << indmap[ndx] << setw(4) << pos
          << setw(5) << L[pos] << endl;
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
    Rcout << "..........k = " << k << "jk = "<< jk << "sz=" << sz << endl;
    int ik = jk;
    for (int i=sz-1;i>=0;i--)
    {
      t[i] = t[i] - L[ik]*L[jk];
      ik++;
    }
  }

  // write results dense matrix t back to L_j
  for (int i=0;i<sz;i++)
  {
     L[posL[i]] = t[i];
  }
}

// just a very simple test how to update L
void mult2(NumericVector& L, double alpha)
{
  const int N = L.size();
  for (int i=0;i<N;i++) {
    L[i] *= alpha;
  }
}

// [[Rcpp::export]]
NumericVector PrintADchol(SEXP arg, NumericVector lambda)
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
  //mult2(L2, 2.0);
  // just print header....
  Rcout << "Entries:" << endl;
  for (int i=0;i<sz;i++) {
    Rcout << setw(3) << i << setw(10) << L[i] << endl;
  }
  Rcout << endl;

  // make set S_j for each column j, see Ng and Peyton:
  vector<set<int> > S = makeSetS(supernodes, rowpointers, colpointers, rowindices);
  const int N = colpointers.size() - 1;
  const int Nsupernodes = supernodes.size()-1;
  for (int J=0; J<Nsupernodes;J++) {
    Rcout << "for supernode " << J << endl;
    map<int, int> indmap = makeIntMap(J, N, rowpointers, rowindices);
    for (int j=supernodes[J];j<supernodes[J+1];j++)
    {
      Rcout << "  for column " << j << endl;
      for (set<int>::const_iterator it=S[j].begin(); it!=S[j].end(); it++)
      {
        int K = *it;
        Cmod_outside_J(L, j, K, indmap,
                       supernodes, rowpointers, colpointers, rowindices);
      }
      Cmod_inside_J(L, j, J, supernodes, colpointers);
      Rcout << "    pivot; cdiv" << endl;
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
  }

  /*
  double sum = 0;
  for (int k=0;k<N;k++)
  {
    int s = colpointers[k];
    sum += 2.0*log(L[s]);
  }
  return sum;
  */
  //Rcout << "check output" << endl;
  //for (int i=0;i<sz;i++)
  //  Rcout << setw(3) << i << setw(12) << L[i] << endl;

  return L;
}


