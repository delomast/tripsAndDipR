#' function to fit beta-binomial model (with or without noise) via EM
#'
#'
#' @keywords internal
#' @noRd

BBpolyEM <- function(refCounts, altCounts, ploidy, h, eps, noise, mdiff, maxrep, maxSubIter){

	lastLLH <- -10000
	repNum <- 0

	nComponents <- ploidy + 1
	if(noise) nComponents <- nComponents + 1
	mix <- rep(1/(nComponents), nComponents)

	tau <- rep(.01, ploidy + 1)

	while(repNum < maxrep){

		# E-step
		mix <- BB_Estep(refCounts, altCounts,
				tau, mix, ploidy, h, eps, noise)

		#  M-step
		optimRes <- stats::optim(
			par = tau,
			fn = llh_calc_BB_Mstep,
			gr = grad_BB_BBnoise,
			method = "L-BFGS-B",
			lower = 10^-7,
			upper = 1 - 10^-7,
			control = list(
				maxit = maxSubIter,
				fnscale = -1
			),
			refCounts = refCounts,
			altCounts = altCounts,
			mixWeights = mix,
			ploidy = ploidy,
			h = h,
			eps = eps
		)
		tau <- optimRes$par

		repNum <- repNum + 1
		# check llh
		if(abs(lastLLH - optimRes$value) < mdiff && repNum > 1) break
		lastLLH <- optimRes$value

	}


	return(list(
		llh = optimRes$value,
		mixProps = mix,
		tau = tau,
		reps = repNum,
		optimRes = optimRes[c("counts", "convergence", "message")]
	))



}
