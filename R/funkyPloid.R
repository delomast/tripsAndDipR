#' calculate LLR's for a group of samples and a given set of ploidy values
#'
#'
#'
#'
#'
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib tripsAndDipR, .registration=TRUE

funkyPloid <- function(counts, counts_alt = NULL, ploidy, h = NULL, eps = NULL,
				   maxRep = 10000, maxDiff = .001){
	# default values
	if(is.null(counts_alt)) {
		counts_alt <- counts[,seq(2, ncol(counts), 2), drop = FALSE]
		counts <- counts[,seq(1,ncol(counts) - 1, 2), drop = FALSE]
		rownames(counts_alt) <- rownames(counts)
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
					mrep = maxRep, mdiff = maxDiff, returnAll = FALSE)
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
