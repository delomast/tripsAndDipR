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

/*
// calculate binomial LH
double binomLH(const double& p, const int& ref, const int& alt){
	return pow(p, ref) * pow(1-p, alt) * exp(logChoose(ref+alt, ref));
}
*/

// calculate binomial log-LH
double binom_log_LH(const double& p, const int& ref, const int& alt){
	return (ref * log(p)) + (alt * log(1-p)) + logChoose(ref+alt, ref);
}

// pmf of multinomial distribution
// @param k a vector of the number of observations in each category
// @param a a vector of the parameters (alpha) describing relative probabilities of each group
double logMultPMF(const vector <double>& k, const vector <double>& a){
	if (k.size() != a.size()) Rcpp::stop("internal error: k not equal to a in logMultPMF");

	int numA = a.size();
	double n = 0;
	for(int i = 0; i < numA; i++) n += k[i];
	double l_sumA = 0;
	for(int i = 0; i < numA; i++) l_sumA += a[i];
	l_sumA = log(l_sumA);

	//start density calculation
	double dens = lgamma(n + 1);
	for(int i = 0; i < numA; i++){
		dens += (k[i] * (log(a[i]) - l_sumA)) - lgamma(k[i] + 1);
	}
	return dens;
}

