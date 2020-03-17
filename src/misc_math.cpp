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

// pmf of Dirichlet-multinomial distribution
// @param k a vector of the number of observations in each category
// @param a a vector of the parameters (alpha) describing the Dirichlet
double logDirichMultPMF(const vector <double>& k, const vector <double>& a){
	if (k.size() != a.size()) Rcpp::stop("internal error: k not equal to a in logdirichMultPMF");
	int numAlpha = a.size();
	double n = 0;
	for(int i = 0; i < numAlpha; i++) n += k[i];
	double sumAlpha = 0;
	for(int i = 0; i < numAlpha; i++) sumAlpha += a[i];

	//start density calculation
	double dens = lgamma(n + 1) + lgamma(sumAlpha) - lgamma(n + sumAlpha);
	for(int i = 0; i < numAlpha; i++){
		dens += lgamma(k[i] + a[i]) - lgamma(k[i] + 1) - lgamma(a[i]);
	}
	return dens;
}

