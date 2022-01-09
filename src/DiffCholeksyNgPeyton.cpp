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
  int nelem = Dim[1];
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
  L["supernodes"] = supernodes;
  L["colpointers"] = rowpointers;
  L["rowpointers"] = colpointers;
  L["rowindices"] =  colindices;
  L["P"] = P_matrix;
  return L;
}


// [[Rcpp::export]]
double PrintADchol(SEXP arg, double lambda)
{
  Rcpp::S4 obj(arg);
  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  NumericMatrix P = Rcpp::clone<Rcpp::NumericMatrix>(obj.slot("P"));

  Rcout << "superNodes:  " <<  supernodes << endl;
  Rcout << "rowpointers: " <<  rowpointers << endl;
  Rcout << "colpointers: " <<  colpointers << endl;
  Rcout << "rowindices:  " <<  rowindices << endl;
  Rcout << endl;
  Rcout << "Explain the structure of sparse matrix:" << endl << endl;
  // in the code below we use variables:
  // s for supernode
  // j for column
  // r for row
  // i for indices
  int Nsupernodes = supernodes.size()-1;
  // for each supernode:
  for (int s=0; s<Nsupernodes;s++) {
    // print current supernode and size:
    Rcout << "SuperNode " << s << " of size "
          << supernodes[s+1] - supernodes[s] << ":" << endl;

    int head = rowpointers[s];
    int end = rowpointers[s+1];
    // for column j:
    for (int j=supernodes[s];j<supernodes[s+1];j++)
    {
      Rcout << " column " << j << ", rows ";
      // get the rows r:
      for (int r=head;r<end;r++)
        Rcout << setw(3) << rowindices[r];
      // get index
      int sc = colpointers[j];
      int ec = colpointers[j+1];
      Rcout << " at indices in entries";
      for (int i=sc;i<ec;i++)
        Rcout << setw(3) << i;
      Rcout << endl;
      // increase head:
      head++;
    }
  }

  return lambda;
}


