#include <Rcpp.h>
#include <math.h>
#include <vector>
#include "misc_math.h"

using namespace std;

// calculate likelihood in M-step of EM-algorithm
// [[Rcpp::export]]
Rcpp::NumericVector llh_calc_BB_Mstep(Rcpp::NumericVector tau, Rcpp::NumericVector refCounts,
                                      Rcpp::NumericVector altCounts,
                                      Rcpp::NumericVector mixWeights,
              int ploidy, Rcpp::NumericVector h, Rcpp::NumericVector eps) {

	// input error checking happens during E-step, so skipping here
	// might be a little dangerous...

	double ploidyD = ploidy;
	unsigned int noisePos = ploidy + 1;

	double nLoci;
	nLoci = refCounts.length();

	bool noise = mixWeights.length() > tau.length();

	/*
	 * calculation of log-likelihood
	 */

	double llh = 0;
	vector <double> oneLocusComps (mixWeights.length());

	for(int i = 0; i < nLoci; i++){

		for(int j=0, max2=ploidyD+1; j < max2; j++ ){
			double scale = (1 - tau[j]) / tau[j];
			double p = j / ploidyD;
			p = (p)*(1 - eps[i]) + (1 - p)*eps[i];
			p = p / ((h[i] * (1 - p)) + p);

			oneLocusComps[j] = logBetaBinomPMF(refCounts[i] + altCounts[i],
                                      refCounts[i], p * scale, (1 - p) * scale) +
                                      	log(mixWeights[j]); // beta-binomials
		}
		if (noise) oneLocusComps[noisePos] = -log(1 + refCounts[i] + altCounts[i]) + log(mixWeights[noisePos]); // uniform noise
		llh += logSumExp(oneLocusComps);
	}

	Rcpp::NumericVector out;
	out.push_back(llh);
	return out;
}
