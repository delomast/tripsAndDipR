#include <Rcpp.h>
#include <math.h>
#include <vector>
#include "misc_math.h"

using namespace std;


// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]


// This runs EM for a given individual and ploidy, then returns the log-likelihood at the maximum
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector genoEM(Rcpp::NumericVector refCounts, Rcpp::NumericVector altCounts,
              int ploidy, Rcpp::NumericVector h, Rcpp::NumericVector eps, int mrep,
              double mdiff, bool returnAll) {

	if (mrep < 1) Rcpp::stop("mrep must be 1 or greater.");
	if (ploidy < 1) Rcpp::stop("ploidy must be 1 or greater.");
	if (refCounts.length() != altCounts.length()) Rcpp::stop("ref and alt counts must have equal lengths.");

	double ploidyD = ploidy; // switch to double for division
	vector <double> gFreq (ploidyD + 1, 1.0/(ploidyD + 1)); // genome-wide freqs
	double nLoci;
	nLoci = refCounts.length();
	// calculations of likelihood for each allele dosage and locus here, using locus specific
	//  eps and h values
	// first index locus, second allele dosage
	// need to add allele specific eps and h values
	vector <vector <double> > pAll;
	for(int i = 0; i < nLoci; i++){
		vector <double> tempVec;
		for(int j=0, max2=ploidyD+1; j < max2; j++ ){
			double p = j / ploidyD;
			p = (p)*(1 - eps[i]) + (1 - p)*eps[i];
			p = p / ((h[i] * (1 - p)) + p);
			tempVec.push_back(binom_log_LH(p, refCounts[i], altCounts[i]));
		}
		pAll.push_back(tempVec);
	}

	// EM to infer genotype category proportions and calculate likelihood
	double llh;
	double lastLlh = -10000; // just preventing division by 0 in first rep
	int repNum = 0;
	for(int rep = 0; rep < mrep; rep++){
		// E- step
		llh = 0;
		vector <double> catCounts (ploidyD + 1, 0);
		for(int i=0; i < nLoci; i++){
			vector <double> genoProbs (ploidyD + 1, 0);
			double rSum = 0;
			for(int j = 0, max2 = ploidyD + 1; j < max2; j++){
				if(gFreq[j] == 0) continue; // prevent log of zero error (genoProbs initiated with 0s)
				genoProbs[j] = exp(pAll[i][j] + log(gFreq[j]));
				rSum += genoProbs[j];
			}
			llh += log(rSum);
			for(int j = 0, max2 = ploidyD + 1; j < max2; j++) catCounts[j] += genoProbs[j] / rSum;
		}
		if ((lastLlh - llh) / lastLlh < mdiff && rep > 0) break;
		lastLlh = llh;

		// M-step
		for(int j = 0, max2 = ploidyD + 1; j < max2; j++) gFreq[j] = catCounts[j] / nLoci;
		repNum++;
	}

	Rcpp::NumericVector returnValues;
	returnValues.push_back(llh);
	if(returnAll){
		// llh, gFreqs, number of reps
		for(int j = 0, max2=ploidyD + 1; j < max2; j++) returnValues.push_back(gFreq[j]);
		returnValues.push_back(repNum);
	}
	return returnValues;
}
