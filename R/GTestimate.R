#' GTestimate
#'
#' @description Normalization of scRNA-seq data using the Simple Good-Turing estimator for relative gene expression estimation.
#' @param object An object containing a scRNA-seq count-matrix in a gene x counts format.
#' @param size.factor Either a numeric vector with one entry per cell, which gives the size factors for cell-level normalization (size factors can be calculated with e.g. scuttle::pooledSizeFactors()) or a number giving the library-size all cells should be scaled to. This defaults to 10,000 in order to behave similar to Seurat's NormalizeData().
#' @param rescale  If TRUE Good-Turing results will be rescaled to add up to 1, before size factors are applied. 
#' @param log1p.transform If TRUE, the GT-estimates are log1p transformed, defaults to TRUE to behave similar to Seurat's NormalizeData()
#' @param assay For Seurat objects: set the assay from which to extract the count matrix, defaults to 'RNA' as this is typically were scRNA-seq count data is stored.
#'
#' For SingleCellExperiment objects: set the assay from which to extract the count matrix, defaults to 'counts' as this is typically were scRNA-seq count data is stored.
#' @param block.size For DelayedArray ojects defines the number of cells read into memory at the same time, this defaults to 1000. Smaller block.sizes require less RAM, but will significantly increase runtime.
#' A minimum of 100 is enforced to protect the user.
#' @details GTestimate is the main function of the GTestimate package.
#' It provides methods for scRNA-seq normalization based on the Good-Turing frequency estimates.
#' GTestimate currently provides methods for regular matrices, sparse dgCMatrix matrices, delayedMatrix matrices, Seurat objects and SingleCellExperiment objects.
#'
#' For matrix input a matrix (either sparse, delayed or dense depending on input) will be returned, containing the normalized gene expression values of the given count-matrix.
#'
#' For Seurat objects a new GTestimate assay will be created within the object (and set as the DefaultAssay). This assay copies the previous counts-data into it's counts slot, and adds the normalized gene expression in the data slot, additionally the missing_mass will be calculated and added as meta-data.
#'
#' For SingleCellExperiment objects a new assay (called GTestimate) containing the normalized gene expression values will be created within the object, additionally the missing_mass will be calculated and added as colData.
#'
#' @return Returns object with normalized gene expression estimates. For Seurat and SingleCellExperiment objects, also adds a Metadata vector missing_mass containing cell-wise estimates of the probability that a next hypothetical unique read would be of a currently unobserved gene.
#' @rdname GTestimate
#' @export
#' @examples
#' library(Seurat)
#' data('pbmc_small')
#' GTestimate(pbmc_small)
GTestimate <- function(object, size.factor = 10000, log1p.transform = TRUE, rescale = FALSE, block.size, assay){
  # Definition of the generic GTestimate method
  if(!(is.numeric(size.factor))) rlang::abort(message = 'Argument size.factor has to be numeric')
  if(!(is.logical(log1p.transform))) rlang::abort(message = 'Argument log1p.transform has to be logical')
  UseMethod("GTestimate", object)
}

#' @method GTestimate matrix
#' @rdname GTestimate
#' @export
GTestimate.matrix <- function(object, size.factor = 10000, log1p.transform = TRUE, rescale = FALSE){
  # Definition of the GTestimate method for dense matrices

  # Detect if method call came from Seurat or SingleCellExperiment object to decide if missing_mass should be returned
  return_mm <- any(c(grepl('GTestimate.Seurat', deparse(sys.calls())),
                     grepl('GTestimate.SingleCellExperiment', deparse(sys.calls()))))

  if(any(object%%1 != 0)) rlang::abort(message = 'Count matrix may only contain integers')
  res_mat <- matrix(0, nrow = nrow(object), ncol = ncol(object), dimnames = dimnames(object))
  for(i in 1:ncol(object)){
    counts <- as.integer(object[,i])
    GT_res <- goodTuring(r = seq.int(from=1L,to=max(counts)), n = tabulate(counts))
    zero <- counts == 0
    m <- match(counts[!zero], names(GT_res))
    res_mat[!zero,i] <- GT_res[m]
  }
  if(return_mm) {
    missing_mass <- 1 - sparseMatrixStats::colSums2(res_mat)
  }
  
  if(rescale){
    res_mat <- t(t(res_mat)/sparseMatrixStats::colSums2(res_mat))
  }
  
  if(length(size.factor) == 1){
    res_mat <- res_mat * size.factor
  } else if(length(size.factor) == ncol(object)){
    res_mat <- t(t(res_mat)*sparseMatrixStats::colSums2(object)/size.factor)
  } else {
    rlang::abort(message = 'Length of size.factor must be 1 or equal to number of columns in object')
  }

  if(log1p.transform) {
    res_mat <- log1p(res_mat)
  }
  if(return_mm){
    return(list(gt_estimates = res_mat, missing_mass = missing_mass))
  } else{
    return(res_mat)
  }

}

#' @method GTestimate dgCMatrix
#' @rdname GTestimate
#' @export
GTestimate.dgCMatrix <- function(object, size.factor = 10000, log1p.transform = TRUE, rescale = FALSE){
  # Definition of the GTestimate method for sparse dgCMatrices, designed to avoid conversion to dense matrix format

  # Detect if method call came from Seurat or SingleCellExperiment object to decide if missing_mass should be returned
  return_mm <- any(c(grepl('GTestimate.Seurat', deparse(sys.calls())),
                     grepl('GTestimate.SingleCellExperiment', deparse(sys.calls()))))

  if(any(object%%1 != 0)) rlang::abort(message = 'Count matrix may only contain integers')
  freq_mat <- sparseMatrixStats::colTabulates(object)[,-1]
  freqs <- as.integer(colnames(freq_mat))
  mat_entries <- object@x
  if (return_mm){
    missing_mass <- c()
    if(length(size.factor) == 1){
      for(i in 1:ncol(object)){
        tmp_range <- seq.int(object@p[i]+1, object@p[i+1])
        GT_res <- goodTuring(r = freqs, n = as.integer(freq_mat[i,]))
        m <- match(object@x[tmp_range], names(GT_res))
        missing_mass <- append(missing_mass, 1 - sum(GT_res[m]))
        if(rescale){
          GT_res[m] <- GT_res[m]/sum(GT_res[m])
        }
        mat_entries[tmp_range]  <- GT_res[m] * size.factor
      }
    } else if(length(size.factor == ncol(object))){
      for(i in 1:ncol(object)){
        tmp_range <- seq.int(object@p[i]+1, object@p[i+1])
        GT_res <- goodTuring(r = freqs, n = as.integer(freq_mat[i,]))
        m <- match(object@x[tmp_range], names(GT_res))
        missing_mass <- append(missing_mass, 1 - sum(GT_res[m]))
        if(rescale){
          GT_res[m] <- GT_res[m]/sum(GT_res[m])
        }
        mat_entries[tmp_range]  <- GT_res[m] * sum(object@x[tmp_range])/ size.factor[i]
      }
    } else {
      rlang::abort(message = paste0('length of size.factor(', length(size.factor),') must be 1 or equal to number of columns in object (',ncol(object), ')'))
    }
  } else{
    if(length(size.factor) == 1){
      for(i in 1:ncol(object)){
        tmp_range <- seq.int(object@p[i]+1, object@p[i+1])
        GT_res <- goodTuring(r = freqs, n = as.integer(freq_mat[i,]))
        m <- match(object@x[tmp_range], names(GT_res))
        if(rescale){
          GT_res[m] <- GT_res[m]/sum(GT_res[m])
        }
        mat_entries[tmp_range]  <- GT_res[m] * size.factor
      }
    } else if(length(size.factor == ncol(object))){
      for(i in 1:ncol(object)){
        tmp_range <- seq.int(object@p[i]+1, object@p[i+1])
        GT_res <- goodTuring(r = freqs, n = as.integer(freq_mat[i,]))
        m <- match(object@x[tmp_range], names(GT_res))
        if(rescale){
          GT_res[m] <- GT_res[m]/sum(GT_res[m])
        }
        mat_entries[tmp_range]  <- GT_res[m] * sum(object@x[tmp_range])/ size.factor[i]
      }
    } else {
      rlang::abort(message = paste0('length of size.factor(', length(size.factor),') must be 1 or equal to number of columns in object (',ncol(object), ')'))
    }
  }
  if (log1p.transform) {
    mat_entries <- log1p(mat_entries)
  } else {
    mat_entries <- mat_entries
  }
  object@x <- mat_entries

  if(return_mm){
    return(list(gt_estimates = object, missing_mass = missing_mass))
  } else{
    return(object)
  }
}

#' @method GTestimate DelayedMatrix
#' @rdname GTestimate
#' @export
GTestimate.DelayedMatrix <- function(object, size.factor = 10000, log1p.transform = TRUE, rescale = FALSE, block.size = 1000L){
  # Definition of the GTestimate method for DelayedMatrices, this processes the input 1000 cells at a time and avoids loading it into memory.

  # Detect if method call came from Seurat or SingleCellExperiment object to decide if missing_mass should be returned
  return_mm <- any(c(grepl('GTestimate.Seurat', deparse(sys.calls())),
                     grepl('GTestimate.SingleCellExperiment', deparse(sys.calls()))))

  if(!(is.numeric(block.size)) | length(block.size) != 1) rlang::abort(message = 'Argument block.size must be a single number')
  grid <- DelayedArray::RegularArrayGrid(refdim = dim(object), spacings = c(nrow(object), min(max(block.size, 100), ncol(object))))
  DelayedArray::setAutoRealizationBackend("HDF5Array")
  sink <- DelayedArray::AutoRealizationSink(dim = dim(object), dimnames = dimnames(object))
  if(return_mm){
    missing_mass <- c()
  }

  if(length(size.factor) == 1){
    for (bid in seq_along(grid)){
      a_viewport <- grid[[bid]]
      block <- read_block(object, a_viewport, as.sparse = NA)
      gt_estimates <- GTestimate::GTestimate(as(block, 'dgCMatrix'), size.factor, log1p.transform, rescale)
      if(return_mm){
        sink <- DelayedArray::write_block(sink = sink, block = gt_estimates$gt_estimates, viewport = a_viewport)
        missing_mass <- c(missing_mass, gt_estimates$missing_mass)
      } else {
        sink <- DelayedArray::write_block(sink = sink, block = gt_estimates, viewport = a_viewport)
      }
    }
  } else if(length(size.factor == ncol(object))){
    sum_of_counts <- 0
    for (bid in seq_along(grid)){
      a_viewport <- grid[[bid]]
      tmp_block <- read_block(object, a_viewport, as.sparse = NA)
      sum_of_counts <- sum_of_counts + sum(tmp_block)
    }
    for (bid in seq_along(grid)){
      a_viewport <- grid[[bid]]
      tmp_range <- start(a_viewport@ranges)[2]:end(a_viewport@ranges)[2]
      block <- read_block(object, a_viewport, as.sparse = NA)
      gt_estimates <- GTestimate::GTestimate(as(block, 'dgCMatrix'), size.factor = size.factor[tmp_range], log1p.transform, rescale)
      if(return_mm){
        sink <- DelayedArray::write_block(sink = sink, block = gt_estimates$gt_estimates, viewport = a_viewport)
        missing_mass <- c(missing_mass, gt_estimates$missing_mass)
      } else {
        sink <- DelayedArray::write_block(sink = sink, block = gt_estimates, viewport = a_viewport)
      }
    }
  } else {
    rlang::abort(message = 'length of size.factor must be equal to number of columns in object')
  }

  DelayedArray::close(sink)
  sink <- as(as(sink, "DelayedArray"), "DelayedMatrix")

  if(return_mm){
    return(list(gt_estimates = sink, missing_mass = missing_mass))
  } else{
    return(sink)
  }
}

#' @method GTestimate Seurat
#' @rdname GTestimate
#' @export
GTestimate.Seurat <- function(object, size.factor = 10000, log1p.transform = TRUE, rescale = FALSE, assay = 'RNA', ...){
  # Definition of the GTestimate method for Seurat objects, extracts count data from the appropriate slot and writes Good-Turing estimates to new assay.
  if(utils::packageVersion('SeuratObject') >= '5.0.0') { # Workaround because Seurat changed slot to layer...
    assay_data <- SeuratObject::GetAssayData(object, layer = 'counts', assay = assay)
  } else {
    assay_data <- SeuratObject::GetAssayData(object, slot = 'counts', assay = assay)
  }

  if (any(dim(assay_data)==c(0,0))){
    rlang::abort(message = 'The chosen assay has no data in the counts slot')
  }
  GT_estimates <- GTestimate(assay_data, size.factor, log1p.transform, rescale, ...)

  object[['GTestimate']] <- SeuratObject::CreateAssayObject(data = GT_estimates$gt_estimates)
  if(utils::packageVersion('SeuratObject') >= '5.0.0') { # Workaround because Seurat changed slot to layer...
    object <- SetAssayData(object, new.data = GetAssayData(object, assay = assay, layer = 'counts'), assay = 'GTestimate', layer = 'counts')
  } else {
    object <- SetAssayData(object, new.data = GetAssayData(object, assay = assay, slot = 'counts'), assay = 'GTestimate', slot = 'counts')
  }

  DefaultAssay(object) <- 'GTestimate'
  object <- SeuratObject::AddMetaData(object, metadata = GT_estimates$missing_mass, col.name = 'missing_mass')
  object <- SeuratObject::LogSeuratCommand(object = object)
  object
}

#' @method GTestimate SingleCellExperiment
#' @rdname GTestimate
#' @export
GTestimate.SingleCellExperiment <- function(object, size.factor = 10000, log1p.transform = TRUE, rescale = FALSE, assay = 'counts', ...){
  # Definition of the GTestimate method for SingleCellExperiment objects, extracts count data from the appropriate slot and writes Good-Turing estimates to new assay.
  assay_data <- SummarizedExperiment::assay(object, assay)
  if (any(dim(assay_data)==c(0,0))){
    rlang::abort(message = 'The chosen assay has no data')
  }
  if (!is.null(sizeFactors(object))) {
    rlang::warn(message = 'The size factors in slot "sizeFactor" of your SingleCellExperiment object were ignored. If you want to use them please set size.factor explicitly e.g. "GTestimate(object, size.factor = sizeFactors(object))". This is done to ensure similar behaviour for Seurat and SingleCellExperiment objects.')
  }
  
  GT_estimates <- GTestimate(assay_data, size.factor, log1p.transform, rescale, ...)
  SummarizedExperiment::assay(object, 'GTestimate') <- GT_estimates$gt_estimates
  SummarizedExperiment::colData(object) <- cbind(SummarizedExperiment::colData(object), missing_mass = GT_estimates$missing_mass)
  object
}
