#' Fit the model for one ploidy
#'
#' For each sample, this returns the log-likelihood and proportions of each genotype category for a given ploidy value
#'
#'
#' @inheritParams funkyPloid
#' @param ploidy The ploidy (as an integer) to fit.
#'
#'
#'
#' @export

genoProps <- function(counts, counts_alt = NULL, ploidy, h = NULL, eps = NULL,
				   maxIter = 10000, maxDiff = .0001, model = c("BB_noise", "BB", "Bin"),
				  maxSubIter = 500){

	model <- match.arg(model)

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
	if(length(ploidy) != 1) stop("Only one ploidy can be specified for this function.")
	if(!isTRUE(all.equal(0, ploidy %% 1))) stop("ploidy must be an integer")

	# fit model
	if(model == "Bin") {
		results <- matrix(NA, nrow = nrow(counts), ncol = ploidy + 3)
	} else if(model == "BB_noise") {
		results <- matrix(NA, nrow = nrow(counts), ncol = (ploidy+2)*2 + 1)
	} else if(model == "BB") {
		results <- matrix(NA, nrow = nrow(counts), ncol = (ploidy+1)*2 + 2)
	}

	mLoci <- rep(NA, nrow(counts))
	for(i in 1:nrow(counts)){
		countBool <- counts[i,] + counts_alt[i,] > 0 # only use loci with reads
		mLoci[i] <- sum(countBool)
		if(mLoci[i] > 0){
			if(model == "Bin"){
				results[i,] <- genoEM(refCounts = counts[i,countBool], altCounts = counts_alt[i,countBool],
					ploidy = ploidy, h = h[countBool], eps = eps[countBool],
					mrep = maxIter, mdiff = maxDiff, returnAll = TRUE)
			} else {
				if (model == "BB_noise") {
					noise <- TRUE
				} else if (model == "BB"){
					noise <- FALSE
				}

				tempRes <- BBpolyEM(refCounts = counts[i,countBool], altCounts = counts_alt[i,countBool],
							ploidy = ploidy, h = h[countBool], eps = eps[countBool],
							noise = noise, mdiff = maxDiff, maxrep = maxIter, maxSubIter = maxSubIter)
				results[i,] <- c(tempRes$llh, tempRes$mixProps, tempRes$tau, tempRes$reps)
			}
		}
	}
	output <- data.frame(Ind = rep(NA, nrow(counts)),
						Loci = mLoci,
					 	numIter = results[,ncol(results)],
						results[,1:(ncol(results)-1), drop = FALSE], stringsAsFactors = FALSE)
	if(!is.null(rownames(counts))){
		output$Ind <- rownames(counts)
	} else {
		output$Ind <- paste0("Row_", 1:nrow(counts))
	}
	colnames(output)[4] <- "LLH"
	colnames(output)[5:(5+ploidy)] <- paste0("ref_", 0:ploidy)
	if(model %in% c("BB_noise", "BB")){
		if(noise){
			colnames(output)[(5+ploidy+1)] <- "noise"
			colnames(output)[(5+ploidy+2):ncol(output)] <- paste0("tau_", 0:ploidy)
		} else {
			colnames(output)[(5+ploidy+1):ncol(output)] <- paste0("tau_", 0:ploidy)
		}
	}

	return(output)
}
