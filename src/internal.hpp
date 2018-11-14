#ifndef INTERNAL_HPP
#define INTERNAL_HPP

#include <Rcpp.h>
#include <cmath>
#include <vector>

double blnratio(double r, double x, double size, double mean, double sd);
Rcpp::NumericVector linspace(const double &start, const double &end, const int &num);
double logistic(const double &x);
double logit(const double &x);
Rcpp::NumericVector logit(const Rcpp::NumericVector &x);
double fxpdf(const double &ratio, const double &x, const double &xc, const double &mean, const double &variance, const double &z);
Rcpp::NumericVector fxpdf(const Rcpp::NumericVector &ratio, const double &x, const double &xc, const double &mean, const double &variance, const double &z);

#endif
