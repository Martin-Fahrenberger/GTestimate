#' GTestimate
#'
#'	Wrapper functions for Good Turing estimation of scRNA-seq relative gene expression estimation.
#' @importFrom rlang abort
#' @importFrom matrixStats colTabulates
#' @export
#' @examples
#' library(Seurat)
#' data('pbmc_small')
#' GTestimate(pbmc_small)


GTestimate <- function(object, ...){
  UseMethod("GTestimate", object)
}
#' @export
GTestimate.matrix <- function(object, ...){
  if(any(object%%1 != 0)) abort(message = 'Count matrix may only contain integers')
  res_mat <- matrix(0, nrow = nrow(object), ncol = ncol(object), dimnames = dimnames(object))
  for(i in 1:ncol(object)){
    counts <- as.integer(object[,i])
    GT_res <- GTestimate::goodTuring(tabulate(counts), seq.int(from=1L,to=max(counts)))
    zero <- counts == 0
    m <- match(counts[!zero], names(GT_res))
    res_mat[!zero,i] <- GT_res[m]
  }
  res_mat
}
#' @export
GTestimate.dgCMatrix <- function(object, ...){
  if(any(object%%1 != 0)) abort(message = 'Count matrix may only contain integers')
  freq_mat <- sparseMatrixStats::colTabulates(object)[,-1]
  freqs <- as.integer(colnames(freq_mat))
  for(i in 1:ncol(object)){
    tmp_range <- seq.int(object@p[i]+1, object@p[i+1])
    GT_res <- GTestimate::goodTuring(as.integer(freq_mat[i,]), freqs)
    m <- match(object@x[tmp_range], names(GT_res))
    #m <- match(counts[counts!=0], names(GT_res))
    #res_mat[counts!=0,i] <- GT_res[m]
    object@x[tmp_range]  <- GT_res[m]
  }
  object
}
#' @export
GTestimate.Seurat <- function(object, ...){
  assay_data <- GetAssayData(object, slot = 'counts')
  GT_estimates <- goodTuringProportions(as.matrix(assay_data))
  object[['GTestimate']] <- CreateAssayObject(data = GT_estimates)
  DefaultAssay(object) <- 'GTestimate'
  object
}
#' @export
GTestimate.SingleCellExperiment <- function(object, ...){
  assay_data <- SummarizedExperiment::assay(object, 'counts')
  GT_estimates <- goodTuringProportions(as.matrix(assay_data))
  SummarizedExperiment::assay(object, 'GTestimate') <- GT_estimates
  object
}
