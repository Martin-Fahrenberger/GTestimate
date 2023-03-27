
# GTestimate

<!-- badges: start -->
<!-- badges: end -->

Implements the simple Good-Turing estimation for scRNA-seq relative gene expression estimation.

## Installation

You can install the development version of GTestimate like so:

``` r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
if (!requireNamespace("sparseMatrixStats", quietly = TRUE))
  devtools::install_github("const-ae/sparseMatrixStats")

devtools::install_github("Martin-Fahrenberger/GTestimate")
```

## Example

This is a basic example of how to use GTestimate to get Good-Turing frequency estimates of the gene expression in a scRNA-seq data-set using Seurat.

``` r
library(GTestimate)
library(Seurat)
data('pbmc_small')
GTestimate(pbmc_small)
```

## Credit

The core implementation of the Good-Turing estimator in C++ has been adapted from Aaron Lun's implementation for the edgeR package.
His implementation was in turn based on Geoffrey Sampson's C code.
