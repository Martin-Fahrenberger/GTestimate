#' GTestimate
#'
#' Implements the simple Good-Turing estimation for scRNA-seq relative gene expression estimation.
#' @param object An object containing a scRNA-seq count-matrix with each column containing one cell's gene expression vector.
#' @param scale.factor Sets the scale factor for cell-level normalization, defaults to 10,000 in order to behave similar to Seurat's RC and LogNormalize
#' @param log1p.transform If TRUE, the GT-estimates are multiplied by the scale.factor and log1p transformed, defaults to TRUE to behave similar to Seurat's LogNornalize
#' @param assay For Seurat objects: set the assay from which to extract the count matrix, defaults to 'RNA' as this is typically were scRNA-seq count data is stored.
#'
#' For SingleCellExperiment objects: set the assay from which to extract the count matrix, defaults to 'counts' as this is typically were scRNA-seq count data is stored.
#' @details GTestimate is the main function of the GTestimate package.
#' It provides methods to calculate the Good-Turing frequency estimates for common scRNA-seq count matrix formats.
#' GTestimate currently provides methods for regular matrices, sparse dgCMatrix matrices, Seurat objects and SingleCellExperiment objects.
#'
#' For matrix input a matrix (either sparse or dense depending on input) will be returned, containing the Good-Turing estimates for the given count-matrix.
#'
#' For Seurat objects a new assay (called GTestimate) will be created within the object (and set as the DefaultAssay), additionally the missing_mass will be calculated and added as meta-data.
#'
#' For SingleCellExperiment objects a new assay (called GTestimate) will be created within the object, additionally the missing_mass will be calculated and added as colData.
#'
#' @return Returns object with Good-Turing gene expression estimates.
#' @rdname GTestimate
#' @export
#' @examples
#' library(Seurat)
#' data('pbmc_small')
#' GTestimate(pbmc_small)


GTestimate <- function(object, scale.factor, log1p.transform, assay){
  UseMethod("GTestimate", object)
}
#' @method GTestimate matrix
#' @rdname GTestimate
#' @export
GTestimate.matrix <- function(object, scale.factor = 10000, log1p.transform = TRUE){
  if(any(object%%1 != 0)) rlang::abort(message = 'Count matrix may only contain integers')
  res_mat <- matrix(0, nrow = nrow(object), ncol = ncol(object), dimnames = dimnames(object))
  for(i in 1:ncol(object)){
    counts <- as.integer(object[,i])
    GT_res <- goodTuring(r = seq.int(from=1L,to=max(counts)), n = tabulate(counts))
    zero <- counts == 0
    m <- match(counts[!zero], names(GT_res))
    res_mat[!zero,i] <- GT_res[m]
  }
  if (log1p.transform) {
    return(log1p(res_mat * scale.factor))
  } else {
    return(res_mat * scale.factor)
  }
}
#' @method GTestimate dgCMatrix
#' @rdname GTestimate
#' @export
GTestimate.dgCMatrix <- function(object, scale.factor = 10000, log1p.transform = TRUE){
  if(any(object%%1 != 0)) rlang::abort(message = 'Count matrix may only contain integers')
  freq_mat <- sparseMatrixStats::colTabulates(object)[,-1]
  freqs <- as.integer(colnames(freq_mat))
  mat_entries <- object@x
  for(i in 1:ncol(object)){
    tmp_range <- seq.int(object@p[i]+1, object@p[i+1])
    GT_res <- goodTuring(r = freqs, n = as.integer(freq_mat[i,]))
    m <- match(object@x[tmp_range], names(GT_res))
    mat_entries[tmp_range]  <- GT_res[m]
  }
  if (log1p.transform) {
    mat_entries <- log1p(mat_entries * scale.factor)
  } else {
    mat_entries <- mat_entries * scale.factor
  }
  object@x <- mat_entries
  return(object)
}
#' @method GTestimate Seurat
#' @rdname GTestimate
#' @export
GTestimate.Seurat <- function(object, scale.factor = 10000, log1p.transform = TRUE, assay = 'RNA'){
  assay_data <- SeuratObject::GetAssayData(object, slot = 'counts', assay = assay)
  if (any(dim(assay_data)==c(0,0))){
    rlang::abort(message = paste0('The chosen assay ', assay, ' has no data in the counts slot'))
  }
  GT_estimates <- GTestimate(assay_data, scale.factor = 1, log1p.transform = F)
  if (log1p.transform){
    object[['GTestimate']] <- SeuratObject::CreateAssayObject(data = log1p(GT_estimates * scale.factor))
  } else{
    object[['GTestimate']] <- SeuratObject::CreateAssayObject(data = GT_estimates * scale.factor)
  }
  DefaultAssay(object) <- 'GTestimate'
  object <- SeuratObject::AddMetaData(object, metadata = 1 - sparseMatrixStats::colSums2(GT_estimates), col.name = 'missing_mass')
  object
}
#' @method GTestimate SingleCellExperiment

#' @rdname GTestimate
#' @export
GTestimate.SingleCellExperiment <- function(object, scale.factor = 10000, log1p.transform = TRUE, assay = 'counts'){
  assay_data <- SummarizedExperiment::assay(object, assay)
  if (any(dim(assay_data)==c(0,0))){
    rlang::abort(message = paste0('The chosen assay ', assay, ' has no data'))
  }
  GT_estimates <- GTestimate(assay_data, scale.factor = 1, log1p.transform = F)
  if (log1p.transform){
    SummarizedExperiment::assay(object, 'GTestimate') <- log1p(GT_estimates * scale.factor)
  } else {
    SummarizedExperiment::assay(object, 'GTestimate') <- GT_estimates * scale.factor
  }
  SummarizedExperiment::colData(object) <- cbind(SummarizedExperiment::colData(object), missing_mass = 1 - sparseMatrixStats::colSums2(GT_estimates))
  object
}
