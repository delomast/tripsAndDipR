#include <Rcpp.h>
#include <math.h>
#include <vector>
#include "misc_math.h"

using namespace std;

/*
 * Calculating the gradient of the log-likelihood of the mixture of beta-binomials
 * With or without noise
 * passed to optim for L-BFGS-B optimization
 */

// derivative of log of beta - binomial wrt alpha
double dl_betabinom_dalpha(const double& n, const double& k, const double& alpha, const double& beta){
	return (R::digamma(k + alpha) - R::digamma(alpha + beta + n)) - (R::digamma(alpha) -
		R::digamma(alpha + beta));
}

// derivative of log of beta - binomial wrt beta
double dl_betabinom_dbeta(const double& n, const double& k, const double& alpha, const double& beta){
	return (R::digamma(n - k + beta) - R::digamma(alpha + beta + n)) - (R::digamma(beta) -
		R::digamma(alpha + beta));
}

// calculating gradient of log-likelihood wrt tau
// each distribution has a separate tau
// [[Rcpp::export]]
Rcpp::NumericVector grad_BB_BBnoise(Rcpp::NumericVector tau, Rcpp::NumericVector refCounts,
                                      Rcpp::NumericVector altCounts,
                                      Rcpp::NumericVector mixWeights,
              int ploidy, Rcpp::NumericVector h, Rcpp::NumericVector eps){


	double ploidyD = ploidy;
	int noisePos = ploidy + 1;

	double nLoci;
	nLoci = refCounts.length();

	bool noise = mixWeights.length() > tau.length();

	vector <double> gradient (tau.length(), 0);

	// for each locus
	for(int i = 0; i < nLoci; i++){

		// log-likelihood for one locus
		double log_likelihood;
		double n = refCounts[i] + altCounts[i];

		vector <double> llh_comps (mixWeights.size());
		for(int t = 0; t < noisePos; t++){
			double p = t / ploidyD;
			p = (p)*(1 - eps[i]) + (1 - p)*eps[i];
			p = p / ((h[i] * (1 - p)) + p);
			double scale = (1 - tau[t]) / tau[t];

			llh_comps[t] = logBetaBinomPMF(n, refCounts[i], p * scale, (1 - p) * scale) + log(mixWeights[t]);
		}
		if (noise) llh_comps[noisePos] = -log(1 + n) + log(mixWeights[noisePos]); // uniform noise
		log_likelihood = logSumExp(llh_comps);

		// for each tau
		for(int t = 0; t < noisePos; t++){
			// only consider the value of p that corresponds to current t b/c others are 0
			double p = t / ploidyD;
			p = (p)*(1 - eps[i]) + (1 - p)*eps[i];
			p = p / ((h[i] * (1 - p)) + p);

			// derivatives of alpha and beta wrt tau
			double dalpha_dt = -p / pow(tau[t], 2);
			double dbeta_dt = -(1 - p) / pow(tau[t], 2);

			double scale = (1 - tau[t]) / tau[t];
			double alpha = p * scale;
			double beta = (1 - p) * scale;

			// dLLH/dalpha * dalpha/dtau + dLLH/dbeta*dbeta/dtau
			gradient[t] += exp(llh_comps[t] - log_likelihood) * ((dl_betabinom_dalpha(n, refCounts[i], alpha, beta) * dalpha_dt) +
				(dl_betabinom_dbeta(n, refCounts[i], alpha, beta) * dbeta_dt));
		}
	}

	// return results
	Rcpp::NumericVector out;
	for (int i = 0, m = gradient.size(); i < m; i++) out.push_back(gradient[i]);

	return out;
}
