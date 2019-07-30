#' Uses read counts for biallelic SNPs to determine if a sample is diploid or triploid
#'
#' \code{tripsAndDip} calculates log-likelihood ratios comparing whether a sample is likely
#' diploid or triploid based on the read counts for biallelic SNPs.
#'
#' \code{tripsAndDip} calculates log-likelihood ratios comparing whether a sample is likely
#' diploid or triploid based on the read counts for biallelic SNPs.This function was designed
#' with amplicon sequencing data in mind, but may be useful for other genotyping techniques
#' that also yield read counts for each allele in a given locus. Full details of the calculations
#' can be found in Delomas (2019) Differentiating diploid and triploid individuals using single
#' nucleotide polymorphisms genotyped by amplicon-sequencing. Molecular Ecology Resources.
#'
#' @param counts Either a numeric matrix or a dataframe with each row corresponding to a different sample.
#'   There are two options for formatting the input. Either
#'   the columns correspond to the read counts for each locus, in a two column per locus format:
#'   column 1 is the read counts for locus1ReferenceAllele, column two is the read counts for locus1AlternateAllele2, locus2Reference, locus2Alternate, ...
#'   OR this contains read counts for the reference allele, and \code{counts_alt} contains read counts for the alternate allele
#'   The rownames should be the sample names.
#' @param counts_alt This is a numeric matrix or a dataframe with each row corresponding to a different sample.
#'   The matrix contains counts for the alternate allele, with samples and loci having the same order as in \code{counts}
#'   If this parameter is NA or NULL, \code{counts} is assumed to have both the reference and alternate allele counts.
#' @param h A numeric vector of h values for each locus in the same order that the loci are ordered in counts.
#'   These h values are as defined by Gerard et al. (2018) "Genotyping polyploids from messy sequencing data"
#'   with h expressed as allele2 / allele1 where allele 1 is listed before allele 2 in counts.
#' @param eps A numeric vector of values for the error rate per read for each locus in the same order that the loci are ordered in counts.
#'   These are expressed as proportions, so a rate of 1\% should be given as 0.01.
#' @param min_reads The minimum number of reads to consider a locus.
#' @param min_loci The minimum number of usable loci in a sample to calculate a log-likelihood ratio.
#' @param binom_p_value The alpha value to use when applying a binomial test to determine
#'   whether to include a locus in the calculation.
#' @return a dataframe with column 1 containing sample names, column 2 containing calculated LLRs (larger means more likely triploid)
#'   and column 3 containing the number of loci used to calculate the LLR
#' @importFrom stats binom.test
#' @examples
#' # make up some data
#' triploid_allele1 <- rbinom(60, 75, 2/3)
#' triploid_allele2 <- 75 - triploid_allele1
#' diploid_allele1 <- rbinom(60, 75, 1/2)
#' diploid_allele2 <- 75 - diploid_allele1
#' # interleave allele counts
#' triploid <- c(rbind(triploid_allele1, triploid_allele2))
#' diploid <- c(rbind(diploid_allele1, diploid_allele2))
#'
#' # create counts matrix
#' allele_counts <- matrix(data = c(triploid, diploid), byrow = TRUE, nrow = 2, ncol = 120)
#' rownames(allele_counts) <- c("triploid", "diploid")
#'
#' #create h and eps vectors
#' h_constant <- rep(1, 60)
#' eps_constant <- rep(.01, 60)
#'
#' #run function
#' ploidy <- tripsAndDip(allele_counts, h = h_constant, eps = eps_constant)
#' @export

tripsAndDip <- function(counts, counts_alt = NA, h, eps, min_reads = 30, min_loci = 15, binom_p_value = .05){
	### input error checking
	if(!is.matrix(counts) && !is.data.frame(counts)){
		stop("Error. counts must be either a matrix or a dataframe.")
	}
	if(is.na(counts_alt) || is.null(counts_alt)){
		if((ncol(counts)/2) != length(h)){
			stop("Error. The number of columns of counts is not equal to twice the length of h.")
		}
	} else {
		if(!is.matrix(counts_alt) && !is.data.frame(counts_alt)){
			stop("Error. counts_alt must be either a matrix or a dataframe.")
		}
		if(ncol(counts) != length(h)){
			stop("Error. The number of columns of counts is not equal to the length of h.")
		}
		if(ncol(counts) != ncol(counts_alt)){
			stop("Error. The number of columns of counts is not equal to the number of columns of counts_alt.")
		}
		if(nrow(counts) != nrow(counts_alt)){
			stop("Error. The number of rows of counts is not equal to the number of rows of counts_alt.")
		}

		# combine counts and counts_alt
		countsNew <- matrix(nrow = nrow(counts), ncol = 2*ncol(counts))
		rownames(countsNew) <- rownames(counts)
		for (i in seq(1, (ncol(countsNew) - 1), 2)){
			countsNew[,i] <- counts[,((i+1)/2)]
			countsNew[,(i+1)] <- counts_alt[,((i+1)/2)]
		}
		counts <- countsNew
	}

	if(length(eps) != length(h)){
		stop("Error. The length of h is not equal to the length of eps.")
	}
	if(min_reads < 1){
		stop("Error. min_reads must be 1 or greater.")
	}
	if(min_loci < 1){
		stop("Error. min_loci must be 1 or greater.")
	}
	if(binom_p_value < 0 || binom_p_value > 1){
		stop("Error. binom_p_value must be between 0 and 1.")
	}

	### calculate llr for each sample
	num_counts <- ncol(counts)
	llr <- apply(counts, 1, function(x){
		#separate allele1 and allele2
		count1 <- x[seq(1, (num_counts - 1), 2)]
		count2 <- x[seq(2, num_counts, 2)]
		#caluculate total read count and which is larger
		n <- count1 + count2
		#determine whether to include loci based on read count
		count1 <- count1[n >= min_reads]
		count2 <- count2[n >= min_reads]
		h_corr <- h[n >= min_reads]
		eps_corr <- eps[n >= min_reads]
		n <- n[n >= min_reads]
		# determine if enough loci
		if (length(n) < min_loci){
			return(c(NA, length(n)))
		}

		k <- mapply(max, count1, count2)
		#flip h as necessary
		h_corr[count2 > count1] <- 1/h_corr[count2 > count1]

		#calculate probabilities
		#based on model with error and allelic bias from Gerard et al. 2018
		p_temp <- (0.6666667*(1 - eps_corr) + (0.3333333)*eps_corr)
		prob_trip <- p_temp / (h_corr*(1 - p_temp) + p_temp)
		p_temp <- 0.5 # (0.5*(1 - eps_corr) + (0.5)*eps_corr) simplifies to 0.5
		prob_dip <- p_temp / (h_corr*(1 - p_temp) + p_temp)

		#determine which loci to include using binomial test
		binom_results <- mapply(function(x,y,z){
				return(binom.test(x, y, z, "greater")$p.value)
			}, k, n, prob_trip)
		binom_results <- binom_results > binom_p_value

		h_corr <- h_corr[binom_results]
		eps_corr <- eps_corr[binom_results]
		n <- n[binom_results]
		k <- k[binom_results]
		prob_trip <- prob_trip[binom_results]
		prob_dip <- prob_dip[binom_results]
		# determine if enough loci
		if (length(n) < min_loci){
			return(c(NA, length(n)))
		}

		#calculate log likelihoods
		trip <- (k*log(prob_trip)) + ((n-k)*log(1-prob_trip))
		dip <- (k*log(prob_dip)) + ((n-k)*log(1-prob_dip))
		#sum accros loci
		trip <- sum(trip)
		dip <- sum(dip)

		return(c(trip - dip, length(n)))

	})
	# associate llr's with sample names and return
	return(data.frame(sample_name = rownames(counts),
			 LLR = llr[1,], loci_used=llr[2,], stringsAsFactors = FALSE, row.names = NULL))
}
