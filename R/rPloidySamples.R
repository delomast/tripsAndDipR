#' simulate read counts for random samples
#'
#' for nSamps samples, read counts are simulated by:
#'   1. If genotypePropsAreKnown is \code{FALSE}, draw genotype proportions from a Dirichlet posterior calculated
#'     using genotypeCounts (treated as a multinomial) and a uniform Dirichlet prior. If
#'     genotypePropsAreKnown is \code{TRUE}, treat genotypeCounts as the true genotype proportions.
#'   2. Draw genotypes from categorical distributions using the genotype proportions from step 1.
#'   3. Draw read counts per locus from a Dirichlet-multinomial with \code{alpha} describing the Dirichlet
#'     and \code{reads} total reads
#'   4. Draw reads of reference and alternative alleles as a binomial using the probability calculated from
#'     the allele dosage of the genotypes drawn in step 2, \code{eps}, and \code{h}
#'
#' @param nSamps number of samples to simulate
#' @param reads number of reads per sample
#' @param truePloidy the true ploidy, as an integer, to simulate
#' @param alpha the alpha values of a Dirichlet distribution describing the
#'   proportion of reads for each locus
#' @param eps Either \code{NULL}, or a numeric vector of values for the error rate per read for each locus in the same
#'   order that the loci are ordered in \code{alpha}
#' @param h Either \code{NULL}, or a numeric vector of h values for each locus in the same order that the
#'   loci are ordered in \code{alpha}
#' @param genotypeCounts a matrix with each row representing a locus (in the same order as \code{alpha}),
#'   and each column representing a genotype with column 1 being 0 copies of the reference allele, column 2
#'   being 1 copy of the reference, ..., column \code{truePloidy + 1} being \code{truePloidy} copies of the reference
#' @param genotypePropsAreKnown boolean indicating whether to threat \code{genotypeCounts} as true proportions,
#'   or the observed number of genotypes in a sample
#' @export

rPloidySamples <- function(nSamps, reads, truePloidy, alpha, eps = NULL, h = NULL, genotypeCounts,
				genotypePropsAreKnown = FALSE){

	# error checking
	if(is.null(h)) h <- rep(1, length(alpha))
	if(is.null(eps)) eps <- rep(.01, length(alpha))
	if(length(h) != length(alpha)) stop("h must be the same length as the number of loci.")
	if(length(h) != length(eps)) stop("h and eps must have the same length.")
	if(length(truePloidy) != 1) stop("Only one ploidy can be specified for this function.")
	if(!isTRUE(all.equal(0, truePloidy %% 1))) stop("truePloidy must be an integer")
	if(ncol(genotypeCounts) != truePloidy + 1) stop("genotypeCounts must have truePloidy + 1 columns")
	if(nrow(genotypeCounts) != length(alpha)) stop("genotypeCounts must have one row for each locus")
	if(nSamps < 1) stop("nSamps must be 1 or greater")

	# make Dirichlet posterior with uniform prior
	if(!genotypePropsAreKnown) genotypeCounts <- genotypeCounts + 1

	allRef <- matrix(nrow = nSamps, ncol = length(alpha))
	allAlt <- matrix(nrow = nSamps, ncol = length(alpha))
	for(r in 1:nSamps){


		# draw proportions for each locus and reads per locus
		rpl <- as.vector(stats::rmultinom(1, reads, ranDirich(alphas = alpha)))

		rRef <- rep(NA, length(alpha))
		for(i in 1:length(alpha)){
			# draw genotype
			if(!genotypePropsAreKnown){
				geno <- sample(0:truePloidy, 1, prob = ranDirich(genotypeCounts[i,]))
			} else {
				geno <- sample(0:truePloidy, 1, prob = genotypeCounts[i,])
			}
			# reads of reference
			pTemp <- geno/truePloidy
			pTemp <- (pTemp)*(1 - eps[i]) + (1 - pTemp)*eps[i]
			pTemp <- pTemp / ((h[i] * (1 - pTemp)) + pTemp);
			rRef[i] <- stats::rbinom(1, rpl[i], pTemp)
		}

		# save read counts
		allRef[r,] <- rRef
		allAlt[r,] <- rpl - rRef
	}
	return(list(
		counts = allRef,
		counts_alt = allAlt
	))
}
