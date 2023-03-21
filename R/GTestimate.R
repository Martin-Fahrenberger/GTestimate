#' GTestimate
#'
#'	Wrapper functions for Good Turing estimation of scRNA-seq relative gene expression estimation.
#' @importFrom rlang abort
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
  mode(object) <- 'integer'
  object_list <- lapply(seq_len(ncol(object)), function(x) (object[,x]))
  freq_table <- matrixStats::colTabulates(object)

  GT_res <- apply(freq_table, 1, function(x) GTestimate::goodTuring(as.table(x)))

  make_exp_vec <- function(props_list, expr_list){
    props_list <- unlist(props_list)
    expr_list <- unlist(expr_list)
    zero <- expr_list == 0
    m <- match(expr_list[!zero], names(props_list))
    expr_list[!zero] <- props_list[m]
    expr_list
  }

  mapply(make_exp_vec, GT_res, object_list)
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
