#include <Rcpp.h>
#include <math.h>
#include <vector>

using namespace std;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// calculate log of the binomial coefficient (choose)
double logChoose(const double& n, const double& k){
	return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
}

// calculate binomial log-LH
double binom_log_LH(const double& p, const int& ref, const int& alt){
	return (ref * log(p)) + (alt * log(1-p)) + logChoose(ref+alt, ref);
}

// log - sum - exp function
// for a set of values x, returns
// log(sum(e^x1 + e^x2 + ...))
double logSumExp(const vector <double>& x){
	if(x.size() < 1) Rcpp::stop("Internal error: vector of length 0 to logSumExp.");
	// find max value
	double maxV = x[0];
	for(int i = 1, max = x.size(); i < max; i++) if (x[i] > maxV) maxV = x[i];

	// calculate and return
	double sum = 0;
	for(int i = 0, max = x.size(); i < max; i++) sum += exp(x[i] - maxV);
	return maxV + log(sum);
}

// calculate log of the beta function
double logBeta(const double& a, const double& b) {
	return lgamma(a) + lgamma(b) - lgamma(a+b);
}

// pmf of beta-binomial distribution
// n is total number of trials
// k is number of successes
// a is alpha
// b is beta
double logBetaBinomPMF(const double& n, const double& k, const double& a, const double& b){
	return logChoose(n, k) + logBeta(k + a, n - k + b) - logBeta(a, b);
}
