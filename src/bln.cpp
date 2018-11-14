#include <Rcpp.h>
#include <cmath>
#include <iostream>
// #include <gsl/gsl_integration.h>
#include "internal.hpp"

#define NUM_INTEGRATE 250

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector blnpdf(const NumericVector &x, const NumericVector &size, const NumericVector &mean = 0, const NumericVector &sd = 1) {
  NumericVector px(x.size());
  NumericVector xc = size - x;
  NumericVector variance = pow(sd, 2);
  NumericVector rd = linspace(0, 1, NUM_INTEGRATE + 1);
  double adjust = (rd[1] - rd[0]) / 2;
  for (int i = 0; i < rd.size(); i++) {
    rd[i] += adjust;
  }
  rd.erase(rd.size() - 1);
  for (int i = 0; i < x.size(); i++) {
    if (variance[i] < 1e-3) {
      px[i] = R::dbinom(x[i], size[i], logistic(mean[i]), false);
    } else {
      double z = lgamma(size[i] + 1) - lgamma(x[i] + 1) - lgamma(xc[i] + 1) - (log(2 * M_PI * variance[i]) * 0.5);
      NumericVector f = exp((log(rd) * (x[i] - 1)) + (log((1 - rd)) * (xc[i] - 1)) - (pow((logit(rd) - mean[i]), 2) / (2 * variance[i])) + z);
      px[i] = exp(-log(NUM_INTEGRATE) + log(sum(f)));
    }
  }
  return px;
}

// [[Rcpp::export]]
NumericVector blnpdfx(const NumericVector &x, const NumericVector &size, const NumericVector &mean = 0, const NumericVector &sd = 1, const bool &drop = false) {
  NumericVector px(x.size());
  NumericVector xc = size - x;
  NumericVector variance = pow(sd, 2);
  NumericVector z(x.size(), 0);
  if (drop) {
    z = lgamma(size + 1) - lgamma(x + 1) - lgamma(xc + 1);
  }
  return px;
}
