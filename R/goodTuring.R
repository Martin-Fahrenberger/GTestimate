#' goodTuring
#'
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
#' @importFrom rlang abort
#' @export
#' @examples
#' test_vec <- c(1,1,2,2,0,0,5,2,0,6,3,0,4,2,1,6,0,7,3,3,4,6,1,3,2,4,5,3,0,0,1)
#' goodTuring(test_vec)

goodTuring <- function(freq_table, conf=1.96)
{
  if(!class(freq_table) == 'table') abort(message('freq_table has to be an object of class "table"'))
  freq_table <- freq_table[freq_table != 0]
  if (is.na(freq_table['0'])){
    n0 <- 0
    n <- as.numeric(freq_table)
    r <- as.numeric(names(freq_table))
  } else{
    n0 <- freq_table['0']
    n <- as.numeric(freq_table[-1])
    r <- as.numeric(names(freq_table[-1]))
  }

#	SGT algorithm (no type checking, as that's enforced above)
	out <- .Call('simple_good_turing', PACKAGE = 'GTestimate', r, n, conf)
	names(out) <- c("P0","proportion")
	names(out$proportion) <- r
	out$proportion
}
