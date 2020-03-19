#' generate random dirichlet
#' @param alphas the alpha defining the Dirichlet
#'
ranDirich <- function(alphas){
	gammas <- sapply(alphas, function(a) rgamma(1,a))
	return(gammas / sum(gammas))
}


tempGrad <- function(tau,
				 refCounts,
			altCounts,
			mixWeights,
			ploidy,
			h,
			eps){

	print("numeric")
	print(numDeriv::grad(llh_calc_BB_Mstep, tau, refCounts = refCounts,
			altCounts = altCounts,
			mixWeights = mixWeights,
			ploidy = ploidy,
			h = h,
			eps = eps))
	print("analytical")
	aGrad <- grad_BB_BBnoise(tau, refCounts,
                                      altCounts,
                                      mixWeights,
              ploidy, h, eps)
	print(aGrad)

	return(aGrad)

}

