#' GTestimate
#'
#'	Wrapper functions for Good Turing estimation of scRNA-seq relative gene expression estimation.
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
  goodTuringProportions(object)
}
#' @export
GTestimate.dgCMatrix <- function(object, ...){
  goodTuringProportions(as.matrix(object))
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
