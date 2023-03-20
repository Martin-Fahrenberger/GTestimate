#' Good Turing Proportions
#'
#' Calculatest the Good Truing proportions for a matrix of expression values
#'
#' @param counts
#'
#' @return matrix of GT proportions
#' @export
#'
#' @examples blabla
#'
goodTuringProportions <- function(counts)
  #	Transform counts to approximately normal expression values
  #	Gordon Smyth
  #	15 Dec 2010.  Last modified 5 Jan 2011.
{
  z <- counts
  nlibs <- ncol(counts)
  for (i in 1:nlibs) {
    g <- goodTuring(counts[,i])
    p0 <- g$P0/g$n0
    zero <- z[,i]==0
    z[zero,i] <- p0
    m <- match(z[!zero,i],g$count)
    z[!zero,i] <- g$proportion[m]
  }
  z
}
