

BBpolyEM <- function(refCounts, altCounts, ploidy, h, eps, noise, mdiff, maxrep){

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
		#  optim
		optimRes <- optim(
			par = tau,
			fn = llh_calc_BB_Mstep,
			gr = tempGrad,
			method = "L-BFGS-B",
			lower = 10^-10,
			upper = 1 - 10^-10,
			control = list(
				maxit = 500,
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
		if((lastLLH - optimRes$value) / lastLLH < mdiff && repNum > 1) break
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
