#ifndef INTERNAL_CPP
#define INTERNAL_CPP

#include "internal.hpp"
#include <Rcpp.h>
#include <cmath>

double blnratio(double r, double x, double size, double mean, double sd) {
  return pow(r, x - 1) *
    pow(1 - r, size - x -1) *
    exp(-(pow((logit(r) - mean), 2) / (2 * pow(sd, 2))));
}

Rcpp::NumericVector linspace(const double &start, const double &end, const int &num) {
  // Thanks to Akavall on StackOverflow
  // https://stackoverflow.com/questions/27028226/python-linspace-in-c
  Rcpp::NumericVector linspaced;
  if (num == 0) {
    return linspaced;
  } else if (num == 1) {
    linspaced.push_back(start);
    return linspaced;
  }
  double delta = (end - start) / (num - 1);
  for (int i = 0; i < num - 1; ++i) {
    linspaced.push_back(start + delta * i);
  }
  linspaced.push_back(end);
  return linspaced;
}

double logistic(const double &x) {
  return 1.0 / (1.0 + exp(-x));
}

double logit(const double &x) {
  return log((x / (1.0 - x)));
}

Rcpp::NumericVector logit(const Rcpp::NumericVector &x) {
  return Rcpp::log((x / (1.0 - x)));
}

double fxpdf(const double &ratio, const double &x, const double &xc, const double &mean, const double &variance, const double &z) {
  return exp((log(ratio) * (x - 1)) + (log(1 - ratio) * (xc - 1)) - ((pow((logit(ratio) - mean), 2) / (2 * variance))) + z);
}

Rcpp::NumericVector fxpdf(const Rcpp::NumericVector &ratio, const double &x, const double &xc, const double &mean, const double &variance, const double &z) {
  Rcpp::NumericVector f(ratio.size());
  f = Rcpp::exp((Rcpp::log(ratio) * (x - 1)) + (Rcpp::log(1 - ratio) * (xc - 1)) - ((Rcpp::pow((logit(ratio) - mean), 2) / (2 * variance))) + z);
  return f;
}

#endif
