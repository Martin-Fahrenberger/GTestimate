
# GTestimate

<!-- badges: start -->
<!-- badges: end -->
GTestimate is a scRNA-seq normalization method. In contrast to other methods it uses the Simple Good-Turing estimator for the per cell gene expression estimation. 
Thereby GTestimate can account for the unobserved genes and avoid overestimation of the observed genes. At default settings it serves as a drop-in replacement for Seurat's NormalizeData().

## Installation

You can install the development version of GTestimate like so:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("sparseMatrixStats", quietly = TRUE))
  devtools::install_github("const-ae/sparseMatrixStats")

devtools::install_github("Martin-Fahrenberger/GTestimate")
```

## Seurat Example

This is a basic example of how to use GTestimate to normalize scRNA-seq data in a Seurat workflow.

``` r
library(GTestimate)
library(Seurat)
data('pbmc_small')

pbmc_small <- GTestimate(pbmc_small) # Instead of NormalizeData(pbmc_small)
pbmc_small <- FindVariableFeatures(pbmc_small)
pbmc_small <- ScaleData(pbmc_small)
pbmc_small <- RunPCA(pbmc_small)
# and so on
```

## SingleCellExperiment Example

This is a basic example of how to use GTestimate to normalize scRNA-seq data in a SingleCellExperiment workflow using size-factors.

```r
library(GTestimate)
library(Seurat)
library(scran)
data('pbmc_small')

pbmc_sce <- as.SingleCellExperiment(pbmc_small)

pbmc_sce <- computeSumFactors(pbmc_sce)
pbmc_sce <- GTestimate(pbmc_sce, size.factors = sizeFactors(pbmc_sce)) # Instead of logNormCounts(sce)
# and so on
```

## Credit

The core implementation of the Simple Good-Turing estimator in C++ has been adapted from Aaron Lun's implementation for the edgeR R-package.
His implementation was in turn based on Geoffrey Sampson's C code acessible at https://www.grsampson.net/D_SGT.c
