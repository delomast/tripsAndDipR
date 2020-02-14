#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <random>
using namespace std;
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// calculate log of the binomial coefficient (choose)
double logChoose(double n, double k){
	return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
}

// calculate binomial LH = need to add in binomial coefficinet?????
// don't need it here b/c same for all, so add it separately and multiply
// by the number of genotypes
double binomLH(const double& p, const int& ref, const int& alt){
	return pow(p, ref) * pow(1-p, alt) * exp(logChoose(ref+alt, ref));
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


//' @export
// [[Rcpp::export]]
double testPloidy(NumericVector counts, NumericVector ploidies, double eps) {
	// just one ploidy right now
	vector <double> k (2,0);
	vector <double> a (2,0);
	a[0] = .75;
	a[1] = .25;
	double ploidy = ploidies[0];
	double llh = 0;
	for(int i=0, max=counts.length(); i < max; i+=2){
		// marginalize true genotype
		double lh = 0;
		for(int j=0, max2=ploidy+1; j < max2; j++ ){
			k[0] = j;
			k[1] = ploidy - j;
			double p = j / ploidy;
			p = (p)*(1 - eps) + (1 - p)*eps;
			lh += binomLH(p, counts[i], counts[i+1]) * exp(logMultPMF(k, a));
			// Rcout << binomLH(p, counts[i], counts[i+1]) << " " << exp(logMultPMF(k, a)) << "\n";
		}
		llh += log(lh);
	}

	return llh;
}

//' @export
// [[Rcpp::export]]
double genoLLHsum(NumericVector counts, NumericVector ploidies, double eps) {
	// just one ploidy right now
	double ploidy = ploidies[0];
	double llh = 0;
	for(int i=0, max=counts.length(); i < max; i+=2){
		// marginalize true genotype
		double lh = 0;
		for(int j=0, max2=ploidy+1; j < max2; j++ ){
			double p = j / ploidy;
			p = (p)*(1 - eps) + (1 - p)*eps;
			lh += binomLH(p, counts[i], counts[i+1]);
		}
		llh += log(lh);
	}

	return llh;
}



//' @export
// [[Rcpp::export]]
double genoEM(NumericVector counts, NumericVector ploidies, double eps, int trep) {

	// just one ploidy right now
	double ploidy = ploidies[0];
	vector <double> gFreq (ploidy+1, 1/(ploidy+1)); // genome-wide freqs
	vector <int> items;
	for(int j=0, max2=ploidy+1; j < max2; j++ ) items.push_back(j);
	int rep = 0;
	while(true){
		// E- step
		vector <double> catCounts (ploidy+1, 0);
		for(int i=0, max=counts.length(); i < max; i+=2){
			vector <double> genoProbs (ploidy+1, 0);
			double rSum = 0;
			for(int j=0, max2=ploidy+1; j < max2; j++ ){
				double p = j / ploidy;
				p = (p)*(1 - eps) + (1 - p)*eps;
				genoProbs[j] = binomLH(p, counts[i], counts[i+1]) * gFreq[j]; // is normalized in sampleC
				rSum += genoProbs[j];
			}
			for(int j=0, max2=ploidy+1; j < max2; j++ ) catCounts[j] += genoProbs[j] / rSum;
		}
		// for(int j=0, max2=ploidy+1; j < max2; j++ ) Rcout << catCounts[j] << " ";
		Rcout << "\n";
		// M-step
		for(int j=0, max2=ploidy+1; j < max2; j++ ) gFreq[j] = catCounts[j] / (counts.length() / 2.0);
		for(int j=0, max2=ploidy+1; j < max2; j++ ) Rcout << gFreq[j] << " ";
		Rcout << "\n";

		// calc llh
		double llh = 0;
		for(int i=0, max=counts.length(); i < max; i+=2){
			double lh = 0;
			for(int j=0, max2=ploidy+1; j < max2; j++){
				double p = j / ploidy;
				p = (p)*(1 - eps) + (1 - p)*eps;
				lh += binomLH(p, counts[i], counts[i+1]) * gFreq[j];
			}
			llh += log(lh);
		}
		Rcout << llh << "\n";
		rep++;
		if(rep > trep) break;
	}

	return 0;
}

