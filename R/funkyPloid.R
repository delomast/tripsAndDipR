#' calculate LLR's for a group of samples and a given set of ploidy values
#'
#'
#'
#'
#'
#'

funkyPloid <- function(counts, counts_alt = NULL, ploidy, h = NULL, eps = NULL,
				   maxRep = 10000, maxDiff = .001){
	# default values
	if(is.null(h)) h <- rep(1, ncol(counts))
	if(is.null(eps)) eps <- rep(.01, ncol(counts))
	if(is.null(counts_alt)) {
		counts_alt <- counts[,seq(2, ncol(counts), 2), drop = FALSE]
		counts <- counts[,seq(1,ncol(counts) - 1, 2), drop = FALSE]
	}

	# input error checking
	if(!all.equal(0, ploidy %% 1)) stop("Ploidies must all be integers")

	# calculate log-likelihoods
	llh <- matrix(NA, nrow = nrow(counts), ncol = length(ploidy))
	for(i in 1:nrow(counts)){
		for(p in 1:length(ploidy)){
			llh[i,p] <- genoEM(refCounts = counts[i,], altCounts = counts_alt[i,],
								ploidy = ploidy[p], h = h, eps = eps, mrep = maxRep,
								mdiff = maxDiff, returnAll = FALSE)
		}
		# calculate log-likelihood ratios comparing model with max likelihood to all other models
		llh[i,] <- max(llh[i,]) - llh[i,]
	}
	llhOutput <- data.frame(Ind = rep(NA, nrow(counts)),
						llh, stringsAsFactors = FALSE)
	if(!is.null(rownames(counts))){
		llhOutput$Ind <- rownames(counts)
	} else {
		llhOutput$Ind <- paste0("Row_", 1:nrow(counts))
	}
	colnames(llhOutput) <- paste0("LLR_max_", ploidy)

	return(llhOutput)
}
