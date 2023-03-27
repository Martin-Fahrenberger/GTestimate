
# GTestimate

<!-- badges: start -->
<!-- badges: end -->

Implements the simple Good-Turing estimation for scRNA-seq relative gene expression estimation.

## Installation

You can install the development version of GTestimate like so:

``` r
# install.packages("devtools")
# devtools::install_github("const-ae/sparseMatrixStats")
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

