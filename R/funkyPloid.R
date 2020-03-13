#' calculate LLR's for a group of samples and a given set of ploidy values
#'
#'
#' @param counts A numeric matrix with each row corresponding to a different sample.
#'   There are two options for formatting the input. Either
#'   the columns correspond to the read counts for each locus, in a two column per locus format:
#'   column 1 is the read counts for locus1ReferenceAllele, column two is the read counts for locus1AlternateAllele2, locus2Reference, locus2Alternate, ...
#'   OR this contains read counts for the reference allele, and \code{counts_alt} contains read counts for the alternate allele
#'   The rownames should be the sample names.
#' @param counts_alt Either \code{NULL} or a numeric matrix with each row corresponding to a different sample.
#'   The matrix contains counts for the alternate allele, with samples and loci having the same order as in \code{counts}
#'   If this parameter is \code{NULL}, \code{counts} is assumed to have both the reference and alternate allele counts.
#' @param ploidy A numeric vector containing the ploidies (as integers) to test.
#' @param h Either \code{NULL}, or a numeric vector of h values for each locus in the same order that the loci are ordered in counts.
#'   These h values are as defined by Gerard et al. (2018) "Genotyping polyploids from messy sequencing data" Genetics 210:789-807.
#'   with h expressed as alternate / reference. These values can be estimated using the R package "updog". If \code{NULL},
#'   h values of 1 (unbiased) are used for all loci.
#' @param eps Either \code{NULL}, or a numeric vector of values for the error rate per read for each locus in the same
#'   order that the loci are ordered in counts.
#'   These are expressed as proportions, so a rate of 1\% should be given as 0.01. These values can be
#'   estimated using the R package "updog". If \code{NULL}, error rates of .01 are assumed for all loci.
#' @param maxIter The maximum number of iterations of the EM algorithm to run for a given sample
#' @param maxDiff This is the maximum proportional change in log-likelihood from the previous iteration to accept
#'   as convergence and stop the EM algorithm.
#'
#'
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib tripsAndDipR, .registration=TRUE
#' @export

funkyPloid <- function(counts, counts_alt = NULL, ploidy, h = NULL, eps = NULL,
				   maxIter = 10000, maxDiff = .001){
	if(!is.matrix(counts)){
		warning("Coercing counts to a matrix.")
		counts <- as.matrix(counts)
	}

	# default values
	if(is.null(counts_alt)) {
		counts_alt <- counts[,seq(2, ncol(counts), 2), drop = FALSE]
		counts <- counts[,seq(1,ncol(counts) - 1, 2), drop = FALSE]
		rownames(counts_alt) <- rownames(counts)
	}

	if(!is.matrix(counts_alt)){
		warning("Coercing counts_alt to a matrix.")
		counts_alt <- as.matrix(counts_alt)
	}

	# input error checking
	if(is.null(h)) h <- rep(1, ncol(counts))
	if(is.null(eps)) eps <- rep(.01, ncol(counts))
	if(length(h) != ncol(counts)) stop("h must be the same length as the number of loci.")
	if(length(h) != length(eps)) stop("h and eps must have the same length.")
	if(any(rownames(counts) != rownames(counts_alt))) warning("counts and counts_alt have different rownames. ",
		"The results will assume that they are in the same order and will use the rownames for counts.")
	if(!isTRUE(all.equal(rep(0, length(ploidy)), ploidy %% 1))) stop("Ploidies must all be integers")

	# calculate log-likelihoods
	llh <- matrix(NA, nrow = nrow(counts), ncol = length(ploidy))
	mLoci <- rep(NA, nrow(counts))
	for(i in 1:nrow(counts)){
		countBool <- counts[i,] + counts_alt[i,] > 0 # only use loci with reads
		mLoci[i] <- sum(countBool)
		if(mLoci[i] > 0){
			for(p in 1:length(ploidy)){
				llh[i,p] <- genoEM(refCounts = counts[i,countBool], altCounts = counts_alt[i,countBool],
					ploidy = ploidy[p], h = h[countBool], eps = eps[countBool],
					mrep = maxIter, mdiff = maxDiff, returnAll = FALSE)
			}
			# calculate log-likelihood ratios comparing model with max likelihood to all other models
			llh[i,] <- max(llh[i,]) - llh[i,]
		}
	}
	llhOutput <- data.frame(Ind = rep(NA, nrow(counts)),
						Loci = mLoci,
						llh, stringsAsFactors = FALSE)
	if(!is.null(rownames(counts))){
		llhOutput$Ind <- rownames(counts)
	} else {
		llhOutput$Ind <- paste0("Row_", 1:nrow(counts))
	}
	colnames(llhOutput)[3:ncol(llhOutput)] <- paste0("LLR_", ploidy)

	return(llhOutput)
}
