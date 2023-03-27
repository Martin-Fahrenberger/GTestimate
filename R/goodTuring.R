#' goodTuring
#'
#'  Wrapper function for the Simple Good-Turing algorithm for frequency estimation
#'	as described by Gale and Sampson (1995).
#'	The called C++ Code is identical to the implementation by Aaron Lun for the edgeR R package
#' @param r: vector of observed frequencies
#' @param n:  vector of frequencies of observed frequencies
#' @param conf: confidence factor for internal fit as described in Gale and Sampson (1995)
#' @return Good-Turing frequency estimates for the frequencies given in r

goodTuring <- function(r, n, conf=1.96)
{
  not_zero <- n != 0
	out <- .Call('simple_good_turing', PACKAGE = 'GTestimate', r[not_zero], n[not_zero], conf)
	names(out) <- r[not_zero]
	out
}
