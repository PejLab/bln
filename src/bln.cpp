#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector blnpdf(const NumericVector &x, const NumericVector &size, const NumericVector &mean = 0, const NumericVector &sd = 1) {
  NumericVector xc = size - x;
  NumericVector px(x.size());
  return xc;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
# timesTwo(42)
*/
