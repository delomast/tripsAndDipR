#' generate random dirichlet
#' @param alphas the alpha defining the Dirichlet
#' @keywords internal
#' @noRd
#'
#'
ranDirich <- function(alphas){
	gammas <- sapply(alphas, function(a) stats::rgamma(1,a))
	return(gammas / sum(gammas))
}
