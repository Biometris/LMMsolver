#include <Rcpp.h>
#include <set>
#include <vector>
#include "AuxFun.h"

using namespace Rcpp;
using namespace std;

// Transform to C++ Notation indices
void transf2C(IntegerVector& ndx)
{
  for (int i=0;i<ndx.size();i++)
  {
    ndx[i]--;
  }
}

