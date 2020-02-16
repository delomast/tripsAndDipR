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


// This runs EM for a given individual and ploidy, then returns the log-likelihood at the maximum
//' @export
// [[Rcpp::export]]
double genoEM(NumericVector counts, int ploidy, double eps, int mrep, double mdiff) {
	if (mrep < 1) Rcpp::stop("mrep must be 1 or greater.");
	double ploidyD = ploidy; // switch to double for division
	vector <double> gFreq (ploidyD+1, 1.0/(ploidyD+1)); // genome-wide freqs
	double nLoci;
	nLoci = counts.length() / 2.0;
	// put calculations of p for each allele dosage and locus here, using locus specific
	//  eps and h values
	// use a vector of vectors with first index locus, second allele dosage
	// actually, just calculate binomLH and save, then iterate to find gFreq
	double llh;
	double lastLlh = -10000; // just preventing division by 0 in first rep
	for(int rep = 0; rep < mrep; rep++){
		// E- step
		llh = 0;
		vector <double> catCounts (ploidyD+1, 0);
		for(int i=0, max=counts.length(); i < max; i+=2){
			vector <double> genoProbs (ploidyD+1, 0);
			double rSum = 0;
			for(int j=0, max2=ploidyD+1; j < max2; j++ ){
				double p = j / ploidyD;
				p = (p)*(1 - eps) + (1 - p)*eps;
				genoProbs[j] = binomLH(p, counts[i], counts[i+1]) * gFreq[j]; // is normalized in sampleC
				rSum += genoProbs[j];
			}
			llh += log(rSum);
			for(int j=0, max2=ploidyD+1; j < max2; j++ ) catCounts[j] += genoProbs[j] / rSum;
		}
		if ((lastLlh - llh) / lastLlh < mdiff && rep > 0) break;
		lastLlh = llh;

		// M-step
		for(int j=0, max2=ploidyD+1; j < max2; j++ ) gFreq[j] = catCounts[j] / nLoci;
	}

	return llh;
}

