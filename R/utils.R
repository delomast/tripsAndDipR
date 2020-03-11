#' generate random dirichlet
#' @param alphas the alpha defining the Dirichlet
#'
ranDirich <- function(alphas){
	gammas <- sapply(alphas, function(a) rgamma(1,a))
	return(gammas / sum(gammas))
}


