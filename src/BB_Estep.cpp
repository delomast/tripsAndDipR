#include <Rcpp.h>
#include <math.h>
#include <vector>
#include "misc_math.h"

using namespace std;

// fitting mixture with overdispersion
// each distribution has a separate tau
// E - step
// [[Rcpp::export]]
Rcpp::NumericVector BB_Estep(Rcpp::NumericVector refCounts, Rcpp::NumericVector altCounts,
                                      Rcpp::NumericVector tau, Rcpp::NumericVector mixWeights,
              int ploidy, Rcpp::NumericVector h, Rcpp::NumericVector eps, bool noise) {

	if (ploidy < 1) Rcpp::stop("ploidy must be 1 or greater.");
	if (refCounts.length() != altCounts.length()) Rcpp::stop("ref and alt counts must have equal lengths.");

	double ploidyD = ploidy;
	unsigned int noisePos = ploidy + 1;

	if (noise && mixWeights.length() != (noisePos + 1)) Rcpp::stop("mixWeights is wrong size.");
	if (!noise && mixWeights.length() != noisePos) Rcpp::stop("mixWeights is wrong size.");
	if (noise && mixWeights.length() != (tau.length() + 1)) Rcpp::stop("mixWeights or tau are the wrong size.");
	if (!noise && mixWeights.length() != tau.length()) Rcpp::stop("mixWeights or tau are the wrong size.");

	double nLoci;
	nLoci = refCounts.length();

	// number of components in the mixture
	int numComponents = noisePos;
	if(noise) numComponents++;


	vector <double> ab (2);
	vector <double> k (2);
	vector <double> oneLocusComps (numComponents);
	vector <double> componentProps (numComponents, 0);

	for(int i = 0; i < nLoci; i++){
		double llh_sum;

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
		if (noise) oneLocusComps[noisePos] = -log(1 + refCounts[i] + altCounts[i]) + log(mixWeights[noisePos]); // uniform noise
		llh_sum = logSumExp(oneLocusComps);

		// distribute to components
		for(int j = 0; j < numComponents; j++) componentProps[j] += exp(oneLocusComps[j] - llh_sum);
	}

	Rcpp::NumericVector out;
	// express as proportions
	for(int j = 0; j < numComponents; j++) out.push_back(componentProps[j] / nLoci);

	return out;
}

