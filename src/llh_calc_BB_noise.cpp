#include <Rcpp.h>
#include <math.h>
#include <vector>
#include "misc_math.h"

using namespace std;

// calculate likelihood
// [[Rcpp::export]]
Rcpp::NumericVector llh_calc_BB_noise(Rcpp::NumericVector refCounts, Rcpp::NumericVector altCounts,
                                      Rcpp::NumericVector tau, Rcpp::NumericVector mixWeights,
              int ploidy, Rcpp::NumericVector h, Rcpp::NumericVector eps) {

	if (ploidy < 1) Rcpp::stop("ploidy must be 1 or greater.");
	if (refCounts.length() != altCounts.length()) Rcpp::stop("ref and alt counts must have equal lengths.");

	double ploidyD = ploidy;
	unsigned int noisePos = ploidy + 1;

	if (mixWeights.length() != noisePos + 1) Rcpp::stop("mixWeights is wrong size.");
	if (mixWeights.length() != tau.length()) Rcpp::stop("mixWeights is wrong size.");

	double nLoci;
	nLoci = refCounts.length();

	/*
	 * calculation of log-likelihood
	 */

	double llh = 0;
	vector <double> ab (2);
	vector <double> k (2);
	vector <double> oneLocusComps (noisePos + 1);

	for(int i = 0; i < nLoci; i++){
		k[0] = refCounts[i];
		k[1] = altCounts[i];
		for(int j=0, max2=ploidyD+1; j < max2; j++ ){
			double scale = (1 - tau[j]) / tau[j];
			double p = j / ploidyD;
			p = (p)*(1 - eps[i]) + (1 - p)*eps[i];
			p = p / ((h[i] * (1 - p)) + p);
			ab[0] = p * scale;
			ab[1] = (1 - p) * scale;
			oneLocusComps[j] = logDirichMultPMF(k, ab) + log(mixWeights[j]); // beta-binomials
		}
		oneLocusComps[noisePos] = -log(1 + refCounts[i] + altCounts[i]) + log(mixWeights[noisePos]); // uniform noise
		llh += logSumExp(oneLocusComps);
	}

	Rcpp::NumericVector out;
	out.push_back(llh);
	return out;
}
