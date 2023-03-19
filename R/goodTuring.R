#' goodTuring
#'	Simple Good-Turing algorithm for frequency estimation
#'	as described by Gale and Sampson (1995)
#'  Sampson's C code translated to C++ and optimized by Aaron Lun
#'	Has been tested against Sampson's C code from
#'	http://www.grsampson.net/RGoodTur.html
#'	and gives identical results.
#'	Gordon Smyth and Aaron Lun
#'	9 Nov 2010.  Last modified 7 Sep 2012.
#' @param x: input expression vector
#' @param conf: confidence value for model switch
#' @return a named list containing
#' P0: the proportion of missing mass
#' proportion: the estimated proportions for the carious count levels
#' count: the count levels
#' n: the number of occurences for each count level
#' n0: the number of 0zeroes in the input
#' @export
#' @examples
#' test_vec <- c(1,1,2,2,5,2,6,3,4,2,6,7,3,3,4,6,3,2,4,5,3)
#' goodTuring(test_vec)

goodTuring <- function(x, conf=1.96)
{
#	Raw frequencies
	x <- as.integer(x)
	if(max(x) < length(x)) {
		n <- tabulate(x+1L)
		n0 <- n[1]
		n <- n[-1]
		max.x <- length(n)
		r <- seq.int(from=1L,to=max.x)
		r <- r[n>0]
		n <- n[n>0]
	} else {
		r <- sort(unique(x))
		z <- match(x,r)
		n <- tabulate(z)
		if(r[1]==0) {
			n0 <- n[1]
			n <- n[-1]
			r <- r[-1]
		} else {
			n0 <- 0
		}
	}

#	SGT algorithm (no type checking, as that's enforced above)
	out <- .Call('simple_good_turing', PACKAGE = 'GTestimate', r, n, conf)
	names(out) <- c("P0","proportion")

	out$count <- r
	out$n <- n
	out$n0 <- n0
	out
}
