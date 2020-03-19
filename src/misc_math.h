// header file for misc_math functions

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

#ifndef MISC_MATH_H
#define MISC_MATH_H

#include <Rcpp.h>
#include <math.h>
#include <vector>

double logChoose(const double& n, const double& k);
double binom_log_LH(const double& p, const int& ref, const int& alt);
double logSumExp(const std::vector <double>& x);
double logBeta(const double& a, const double& b);
double logBetaBinomPMF(const double& n, const double& k, const double& a, const double& b);

#endif
