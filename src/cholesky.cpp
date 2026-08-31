// sparse left-looking Cholesky using supernodes.
//
// For details see:
// Ng, Esmond G., and Barry W. Peyton.,
// "Block sparse Cholesky algorithms on advanced uniprocessor computers."
// SIAM Journal on Scientific Computing 14, no. 5 (1993): 1034-1056.
//
// Furrer, Reinhard, and Stephan R. Sain.
// "spam: A sparse matrix R package with emphasis on MCMC
// methods for Gaussian Markov random fields."
// Journal of Statistical Software 36 (2010): 1-25.

//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <set>
#include <vector>
#include "AuxFun.h"
#include "SparseMatrix.h"
#include "cholesky.h"

#include <chrono>

using namespace std::chrono;

using namespace Rcpp;
using namespace std;

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

// [[Rcpp::export]]
void arma_info()
{
  Rcpp::Rcout << "Armadillo version: "
              << arma::arma_version::as_string() << "\n";

  Rcpp::Rcout << "BLAS enabled: "
              << arma::arma_config::blas << "\n";

  Rcpp::Rcout << "LAPACK enabled: "
              << arma::arma_config::lapack << "\n";

  Rcpp::Rcout << "BLAS wrapper: "
              << arma::arma_config::wrapper << "\n";
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

// Adjust column j for all columns in supernode K:
void cmod2(NumericVector& L, int j, int K, int sz,
           NumericVector& t,
           const IntegerVector& indmap,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  // init t:
  for (int i=0;i<sz;i++)
  {
    t[i] = 0.0;
  }

  const int& sCol = supernodes[K];
  const int& eCol = supernodes[K+1];

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

// [[Rcpp::export]]
void chol_last_supernode(
    NumericVector& L,
    const IntegerVector& supernodes,
    const IntegerVector& colpointers)
{
  const int J = supernodes.size()-2;
  //Rcout << "J = " << J << endl;

  const int j0 = supernodes[J];
  const int j1 = supernodes[J + 1];

  const int n = j1 - j0;       // should be 492

  // ------------------------------------------------------------
  // 1. Sparse -> dense
  // ------------------------------------------------------------
  auto s0 = high_resolution_clock::now();

  arma::mat A(n, n, arma::fill::zeros);

  double* l = L.begin();

  for (int j = 0; j < n; ++j)
  {
    const int col = j0 + j;
    const int s = colpointers[col];

    for (int i = j; i < n; ++i)
      A(i, j) = A(j, i) = l[s + (i - j)];
  }

  auto s1 = high_resolution_clock::now();

  double time_s2d =
    duration<double>(s1 - s0).count();

  // ------------------------------------------------------------
  // 2. Full dense Cholesky
  // ------------------------------------------------------------
  auto s2 = high_resolution_clock::now();

  A = arma::chol(A).t();

  auto s3 = high_resolution_clock::now();

  double time_chol = duration<double>(s3 - s2).count();


  // ------------------------------------------------------------
  // 3. Dense -> sparse
  // ------------------------------------------------------------
  auto s4 = high_resolution_clock::now();

  for (int j = 0; j < n; ++j)
  {
    const int col = j0 + j;
    const int s = colpointers[col];

    for (int i = j; i < n; ++i)
      l[s + (i - j)] = A(i, j);
  }

  auto s5 = high_resolution_clock::now();

  double time_d2s =
    duration<double>(s5 - s4).count();


  // ------------------------------------------------------------
  // Total
  // ------------------------------------------------------------
  double time_total = duration<double>(s5 - s0).count();

  return;
}


void cholesky(NumericVector& L,
           const IntegerVector& supernodes,
           const IntegerVector& rowpointers,
           const IntegerVector& colpointers,
           const IntegerVector& rowindices)
{
  auto start_chol = high_resolution_clock::now();

  const int N = colpointers.size() - 1;
  const int Nsupernodes = supernodes.size()-1;

  // linked lists, see section 4.2 Ng and Peyton
  IntegerVector HEAD(N,-1);
  IntegerVector LINK(Nsupernodes,-1);

  IntegerVector colhead = clone(rowpointers);
  for (int J=0; J<Nsupernodes;J++)
  {
    int szNode = supernodes[J+1] - supernodes[J];
    colhead[J] += szNode-1;
    if (colhead[J] < rowpointers[J+1]-1)
    {
      int rNdx = rowindices[colhead[J]+1];
      insert(HEAD, LINK, rNdx, J);
    }
  }

  IntegerVector indmap(N,0);
  NumericVector t(N);
  int last_sn = Nsupernodes-1;

  // for each supernode J
  for (int J=0; J<Nsupernodes;J++) {


    long long calls = 0;
    long long work = 0;

    int sz_node = supernodes[J+1] - supernodes[J];

    // Phase 1
    makeIndMap(indmap, J, rowpointers, rowindices);
    auto start_cmod2 = high_resolution_clock::now();

    for (int j=supernodes[J];j<supernodes[J+1];j++)
    {
      int K = HEAD[j];

      while (K!=-1)
      {

        int nextK = LINK[K];
        int sz = rowpointers[K+1] - colhead[K];

        const int wK = supernodes[K+1] - supernodes[K];
        calls++;
        work += static_cast<long long>(sz) * wK;
        cmod2(L, j, K, sz, t, indmap, supernodes, rowpointers, colpointers, rowindices);

        colhead[K]++;
        if (colhead[K] < rowpointers[K+1])
        {
          int rNdx = rowindices[colhead[K]];
          insert(HEAD, LINK, rNdx, K);
        }
        K = nextK;
      }
      HEAD[j] = -1;
    }

    auto end_cmod2 = high_resolution_clock::now();
    double elapsed_cmod2 = duration<double>(end_cmod2 - start_cmod2).count();

    // update
    colhead[J]++;

    NumericVector copy_L = clone(L);

    if (J == last_sn) {
      Rcout << "Comparison using arma vs original code" << endl;

      auto start_chol_last = high_resolution_clock::now();

      chol_last_supernode(copy_L, supernodes, colpointers);

      auto end_chol_last = high_resolution_clock::now();

      double elapsed_chol_last = duration<double>(end_chol_last - start_chol_last).count();

      Rcout << "Supernode " << setw(3) << J
            << " sznode: " << setw(3) << sz_node
            << " cmod2:  " << setw(3) << elapsed_cmod2
            << " cmod_arma: " << setw(4) << elapsed_chol_last
            << std::endl;
    }

    // Phase 2
    auto start = high_resolution_clock::now();
    for (int j=supernodes[J];j<supernodes[J+1];j++) {
      cmod1(L, j, J, supernodes, colpointers);
      cdiv(L, j, colpointers);
    }
    auto end = high_resolution_clock::now();
    double elapsed = duration<double>(end - start).count();
    if (J == last_sn) {
      if (!Rcpp::is_true(Rcpp::all(Rcpp::abs(copy_L - L) < 1e-09)))
        Rcpp::stop("copy_L and L differ");

    }

    if (sz_node > 2) {
      Rcout << "Supernode " << setw(3) << J << " sznode: " << setw(3) << sz_node <<
        " cmod2: " << setw(4) << elapsed_cmod2 << " cmod1/cdiv: " << setw(4) << elapsed << endl;
    }

  }
  auto end_chol = high_resolution_clock::now();
  double elapsed_total = duration<double>(end_chol - start_chol).count();
  Rcout << endl << "Total time Cholesky: " << elapsed_total << endl << endl;
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


/*
 *
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

*/

/*
// Just to show the structure of the sparse cholesky matrix with supernodes.
// [[Rcpp::export]]
NumericMatrix PrintCholesky(Rcpp::S4 obj)
{
  Rcout << "Class: " << as<std::string>(obj.attr("class")) << std::endl;

  IntegerVector supernodes = obj.slot("supernodes");
  IntegerVector colpointers = obj.slot("colpointers");
  IntegerVector rowpointers = obj.slot("rowpointers");
  IntegerVector rowindices = obj.slot("rowindices");
  IntegerVector pivot = obj.slot("pivot");
  IntegerVector invpivot = obj.slot("invpivot");

  NumericVector L = obj.slot("entries");

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
