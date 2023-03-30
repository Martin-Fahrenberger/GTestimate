#' GTestimate
#'
#' Implements the simple Good-Turing estimation for scRNA-seq relative gene expression estimation.
#' @param object An object containing a scRNA-seq count-matrix with each column containing one cell's gene expression vector.
#' @param scale.factor Sets the scale factor for cell-level normalization, defaults to 10,000 in order to behave similar to Seurat's RC and LogNormalize
#' @param log1p.transform If TRUE, the GT-estimates are multiplied by the scale.factor and log1p transformed, defaults to TRUE to behave similar to Seurat's LogNornalize
#' @param assay For Seurat objects: set the assay from which to extract the count matrix, defaults to 'RNA' as this is typically were scRNA-seq count data is stored.
#'
#' For SingleCellExperiment objects: set the assay from which to extract the count matrix, defaults to 'counts' as this is typically were scRNA-seq count data is stored.
#' @param block_size For DelayedArray ojects defines the number of cells read into memory at the same time, this defaults to 1000. Smaller block_sizes require less RAM, but will significantly increase runtime.
#' A minimum of 100 is enforced to protect the user.
#' @details GTestimate is the main function of the GTestimate package.
#' It provides methods to calculate the Good-Turing frequency estimates for common scRNA-seq count matrix formats.
#' GTestimate currently provides methods for regular matrices, sparse dgCMatrix matrices, delayedMatrix matrices, Seurat objects and SingleCellExperiment objects.
#'
#' For matrix input a matrix (either sparse, delayed or dense depending on input) will be returned, containing the Good-Turing estimates for the given count-matrix.
#'
#' For Seurat objects a new assay (called GTestimate) will be created within the object (and set as the DefaultAssay), additionally the missing_mass will be calculated and added as meta-data.
#'
#' For SingleCellExperiment objects a new assay (called GTestimate) will be created within the object, additionally the missing_mass will be calculated and added as colData.
#'
#' @return Returns object with Good-Turing gene expression estimates. For Seurat and SingleCellExperiment objects also adds a Metadata vector missing_mass which gives the estimated portion of mRNA's in the cell belonging to unobserved genes (refered to as P0 in Gale and Sampson 1995)
#' @rdname GTestimate
#' @export
#' @examples
#' library(Seurat)
#' data('pbmc_small')
#' GTestimate(pbmc_small)
GTestimate <- function(object, scale.factor, log1p.transform, assay, block_size){
  # Definition of the generic GTestimate method
  UseMethod("GTestimate", object)
}

#' @method GTestimate matrix
#' @rdname GTestimate
#' @export
GTestimate.matrix <- function(object, scale.factor = 10000, log1p.transform = TRUE){
  # Definition of the GTestimate method for dense matrices
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
  # Definition of the GTestimate method for sparse dgCMatrices, designed to avoid conversion to dense matrix format
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

#' @method GTestimate DelayedMatrix
#' @rdname GTestimate
#' @export
GTestimate.DelayedMatrix <- function(object, scale.factor = 10000, log1p.transform = TRUE, block_size = 1000L){
  # Definition of the GTestimate method for DelayedMatrices, this processes the input 1000 cells at a time and avoids loading it into memory.
  grid <- DelayedArray::RegularArrayGrid(refdim = dim(object), spacings = c(nrow(object), min(max(block_size, 100), ncol(object))))
  DelayedArray::setAutoRealizationBackend("HDF5Array")
  sink <- DelayedArray::AutoRealizationSink(dim = dim(object), dimnames = dimnames(object))
  for (bid in seq_along(grid)){
    a_viewport <- grid[[bid]]
    block <- read_block(object, a_viewport, as.sparse = NA)
    block <- GTestimate::GTestimate(as(block, 'dgCMatrix'), scale.factor, log1p.transform)
    sink <- DelayedArray::write_block(sink = sink, block = block, viewport = a_viewport)
  }
  DelayedArray::close(sink)
  return(as(as(sink, "DelayedArray"), "DelayedMatrix"))
}

#' @method GTestimate Seurat
#' @rdname GTestimate
#' @export
GTestimate.Seurat <- function(object, scale.factor = 10000, log1p.transform = TRUE, assay = 'RNA', block_size = 1000L){
  # Definition of the GTestimate method for Seurat objects, extracts count data from the appropriate slot and writes Good-Turing estimates to new assay.
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
  object <- LogSeuratCommand(object = object)
  object
}

#' @method GTestimate SingleCellExperiment
#' @rdname GTestimate
#' @export
GTestimate.SingleCellExperiment <- function(object, scale.factor = 10000, log1p.transform = TRUE, assay = 'counts', block_size = 1000L){
  # Definition of the GTestimate method for SingleCellExperiment objects, extracts count data from the appropriate slot and writes Good-Turing estimates to new assay.
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
