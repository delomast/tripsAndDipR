#' calculate proportions of each genotype category for a group of samples and a given ploidy value
#'
#'
#'
#'
#'
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib tripsAndDipR, .registration=TRUE

genoProps <- function(counts, counts_alt = NULL, ploidy, h = NULL, eps = NULL,
				   maxIter = 10000, maxDiff = .001){
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
	if(length(ploidy) != 1) stop("Only one ploidy can be specified for this function.")
	if(!isTRUE(all.equal(0, ploidy %% 1))) stop("Ploidies must be an integer")

	# fit model
	results <- matrix(NA, nrow = nrow(counts), ncol = ploidy + 3)
	mLoci <- rep(NA, nrow(counts))
	for(i in 1:nrow(counts)){
		countBool <- counts[i,] + counts_alt[i,] > 0 # only use loci with reads
		mLoci[i] <- sum(countBool)
		if(mLoci[i] > 0){
			results[i,] <- genoEM(refCounts = counts[i,countBool], altCounts = counts_alt[i,countBool],
				ploidy = ploidy, h = h[countBool], eps = eps[countBool],
				mrep = maxIter, mdiff = maxDiff, returnAll = TRUE)
		}
	}
	results <-
	output <- data.frame(Ind = rep(NA, nrow(counts)),
						Loci = mLoci,
					 	numIter = results[,ncol(results)],
						results[,1:(ncol(results)-1)], stringsAsFactors = FALSE)
	if(!is.null(rownames(counts))){
		output$Ind <- rownames(counts)
	} else {
		output$Ind <- paste0("Row_", 1:nrow(counts))
	}
	colnames(output)[4] <- "LLH"
	colnames(output)[5:ncol(output)] <- paste0("ref_", 0:ploidy)

	return(output)
}
