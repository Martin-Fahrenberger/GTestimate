library(Seurat)
library(SingleCellExperiment)
library(DelayedArray)
library(HDF5Array)
data('pbmc_small')

test_dgCMatrix <- GetAssayData(pbmc_small, slot = 'counts')
test_matrix <- as.matrix(test_dgCMatrix)
test_sf <- scuttle::computePooledFactors(as.SingleCellExperiment(pbmc_small))$sizeFactor

test_delayed <- as(as(test_matrix, 'HDF5Matrix'), 'DelayedMatrix')

test_seurat_dense <- pbmc_small
test_seurat_dense@assays$RNA@counts <- test_matrix

test_sce <- as.SingleCellExperiment(pbmc_small)
test_sce_delayed <- test_sce
assay(test_sce_delayed, 'counts') <- test_delayed

{
  test_that("GTestimate.matrix size.factor = 1 log1p.transform = FALSE", {
    testthat::expect_equal(
      edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0),
      GTestimate::GTestimate(test_matrix, size.factor = 1, log1p.transform = F)
    )
  })

  test_that("GTestimate.matrix size.factor = 10000 log1p.transform = FALSE", {
    testthat::expect_equal(
      (edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0),
      GTestimate::GTestimate(test_matrix, size.factor = 10000, log1p.transform = F)
    )
  })

  test_that("GTestimate.matrix size.factor = 10000 log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p((edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(test_matrix, size.factor = 10000, log1p.transform = T)
    )
  })

  test_that("GTestimate.matrix default parameters", {
    testthat::expect_equal(
      log1p((edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(test_matrix)
    )
  })

  test_that("GTestimate.matrix size.factor = computePooledFactors() log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix)%>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      GTestimate::GTestimate(test_matrix, size.factor = test_sf, log1p.transform = T)
    )
  })
}

{
  test_that("GTestimate.dgCMatrix size.factor = 1 log1p.transform = FALSE", {
    testthat::expect_equal(
      edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0),
      as.matrix(GTestimate::GTestimate(test_dgCMatrix, size.factor = 1, log1p.transform = F))
    )
  })

  test_that("GTestimate.dgCMatrix size.factor = 10000 log1p.transform = FALSE", {
    testthat::expect_equal(
      (edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0),
      as.matrix(GTestimate::GTestimate(test_dgCMatrix, size.factor = 10000, log1p.transform = F))
    )
  })

  test_that("GTestimate.dgCMatrix size.factor = 10000 log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p((edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0)),
      as.matrix(GTestimate::GTestimate(test_dgCMatrix, size.factor = 10000, log1p.transform = T))
    )
  })

  test_that("GTestimate.dgCMatrix size.factor = computePooledFactors() log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(GTestimate::GTestimate(test_dgCMatrix, size.factor = test_sf, log1p.transform = T))
    )
  })
}
{
  test_that("GTestimate.DelayedMatrix size.factor = 1 log1p.transform = FALSE", {
    testthat::expect_equal(
      edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0),
      as.matrix(GTestimate::GTestimate(test_delayed, size.factor = 1, log1p.transform = F))
    )
  })

  test_that("GTestimate.DelayedMatrix size.factor = 10000 log1p.transform = FALSE", {
    testthat::expect_equal(
      (edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0),
      as.matrix(GTestimate::GTestimate(test_delayed, size.factor = 10000, log1p.transform = F))
    )
  })

  test_that("GTestimate.DelayedMatrix size.factor = 10000 log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p((edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0)),
      as.matrix(GTestimate::GTestimate(test_delayed, size.factor = 10000, log1p.transform = T))
    )
  })

  test_that("GTestimate.DelayedMatrix size.factor = computePooledFactors() log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(GTestimate::GTestimate(test_delayed, size.factor = test_sf, log1p.transform = T))
    )
  })
}

{
  test_that("GTestimate.Seurat size.factor = 1 log1p.transform = FALSE", {
    testthat::expect_equal(
      edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0),
      as.matrix(GetAssayData(GTestimate::GTestimate(pbmc_small, size.factor = 1, log1p.transform = F), slot = 'data', assay = 'GTestimate'))
    )
  })

  test_that("GTestimate.Seurat size.factor = 10000 log1p.transform = FALSE", {
    testthat::expect_equal(
      (edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0),
      as.matrix(GetAssayData(GTestimate::GTestimate(pbmc_small, size.factor = 10000, log1p.transform = F), slot = 'data', assay = 'GTestimate'))
    )
  })

  test_that("GTestimate.Seurat size.factor = 10000 log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p((edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0)),
      as.matrix(GetAssayData(GTestimate::GTestimate(pbmc_small, size.factor = 10000, log1p.transform = T), slot = 'data', assay = 'GTestimate'))
    )
  })

  test_that("GTestimate.Seurat size.factor = computePooledFactors() log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(GetAssayData(GTestimate::GTestimate(pbmc_small, size.factor = test_sf, log1p.transform = T), slot = 'data', assay = 'GTestimate'))
    )
  })

  test_that("GTestimate.Seurat size.factor = 1 log1p.transform = FALSE missing mass calculation", {
    testthat::expect_equal(
      1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(pbmc_small, size.factor = 1, log1p.transform = F)$missing_mass
    )
  })

  test_that("GTestimate.Seurat size.factor = 10000 log1p.transform = TRUE missing mass calculation", {
    testthat::expect_equal(
      1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(pbmc_small, size.factor = 10000, log1p.transform = T)$missing_mass
    )
  })

  test_that("GTestimate.Seurat size.factor = computePooledFactors() log1p.transform = TRUE missing mass calculation", {
    testthat::expect_equal(
      1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(pbmc_small, size.factor = test_sf, log1p.transform = T)$missing_mass
    )
  })

  test_that("GTestimate.Seurat dense size.factor = computePooledFactors() log1p.transform = TRUE missing mass calculation", {
    testthat::expect_equal(
      1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(test_seurat_dense, size.factor = test_sf, log1p.transform = T)$missing_mass
    )
  })
}
{
  test_that("GTestimate.SingleCellExperiment size.factor = 1 log1p.transform = FALSE", {
    testthat::expect_equal(
      edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0),
      as.matrix(GTestimate::GTestimate(test_sce, size.factor = 1, log1p.transform = F)@assays@data$GTestimate)
      )
  })

  test_that("GTestimate.SingleCellExperiment size.factor = 10000 log1p.transform = FALSE", {
    testthat::expect_equal(
      (edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0),
      as.matrix(assay(GTestimate::GTestimate(test_sce, size.factor = 10000, log1p.transform = F), 'GTestimate'))
    )
  })

  test_that("GTestimate.SingleCellExperiment size.factor = 10000 log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p((edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0)),
      as.matrix(assay(GTestimate::GTestimate(test_sce, size.factor = 10000, log1p.transform = T), 'GTestimate'))
    )
  })

  test_that("GTestimate.SingleCellExperiment size.factor = computePooledFactors() log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(assay(GTestimate::GTestimate(test_sce, size.factor = test_sf, log1p.transform = T), 'GTestimate'))
    )
  })

  test_that("GTestimate.SingleCellExperiment delayed size.factor = computePooledFactors() log1p.transform = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(assay(GTestimate::GTestimate(test_sce_delayed, size.factor = test_sf, log1p.transform = T), 'GTestimate'))
    )
  })

  test_that("GTestimate.SingleCellExperiment size.factor = 1 log1p.transform = FALSE missing mass calculation", {
    testthat::expect_equal(
      as.numeric(1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      GTestimate::GTestimate(test_sce, size.factor = 1, log1p.transform = F)$missing_mass
    )
  })

  test_that("GTestimate.SingleCellExperiment size.factor = 10000 log1p.transform = TRUE missing mass calculation", {
    testthat::expect_equal(
      as.numeric(1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      GTestimate::GTestimate(test_sce, size.factor = 10000, log1p.transform = T)$missing_mass
    )
  })

  test_that("GTestimate.SingleCellExperiment size.factor = computePooledFactors() log1p.transform = TRUE missing mass calculation", {
    testthat::expect_equal(
      as.numeric(1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      GTestimate::GTestimate(test_sce, size.factor = test_sf, log1p.transform = T)$missing_mass
    )
  })

  test_that("GTestimate.SingleCellExperiment delayed size.factor = computePooledFactors() log1p.transform = TRUE missing mass calculation", {
    testthat::expect_equal(
      as.numeric(1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      GTestimate::GTestimate(test_sce_delayed, size.factor = test_sf, log1p.transform = T)$missing_mass
    )
  })
}

{
  test_that("GTestimate.matrix size.factor = 1 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      GTestimate::GTestimate(test_matrix, size.factor = 1, log1p.transform = F, rescale = T)
    )
  })
  
  test_that("GTestimate.matrix size.factor = 10000 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))) * 10000,
      GTestimate::GTestimate(test_matrix, size.factor = 10000, log1p.transform = F, rescale = T)
    )
  })
  
  test_that("GTestimate.matrix size.factor = 10000 log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))) * 10000),
      GTestimate::GTestimate(test_matrix, size.factor = 10000, log1p.transform = T, rescale = T)
    )
  })
  
  test_that("GTestimate.matrix default parameters", {
    testthat::expect_equal(
      log1p((edgeR::goodTuringProportions(test_matrix) * 10000) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(test_matrix)
    )
  })
  
  test_that("GTestimate.matrix size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix)%>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix)%>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      GTestimate::GTestimate(test_matrix, size.factor = test_sf, log1p.transform = T, rescale = T)
    )
  })
}

{
  test_that("GTestimate.dgCMatrix size.factor = 1 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      as.matrix(GTestimate::GTestimate(test_dgCMatrix, size.factor = 1, log1p.transform = F, rescale = T))
    )
  })
  
  test_that("GTestimate.dgCMatrix size.factor = 10000 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * 10000),
      as.matrix(GTestimate::GTestimate(test_dgCMatrix, size.factor = 10000, log1p.transform = F, rescale = T))
    )
  })
  
  test_that("GTestimate.dgCMatrix size.factor = 10000 log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))) * 10000),
      as.matrix(GTestimate::GTestimate(test_dgCMatrix, size.factor = 10000, log1p.transform = T, rescale = T))
    )
  })
  
  test_that("GTestimate.dgCMatrix size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(GTestimate::GTestimate(test_dgCMatrix, size.factor = test_sf, log1p.transform = T, rescale = T))
    )
  })
}
{
  test_that("GTestimate.DelayedMatrix size.factor = 1 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      as.matrix(GTestimate::GTestimate(test_delayed, size.factor = 1, log1p.transform = F, rescale = T))
    )
  })
  
  test_that("GTestimate.DelayedMatrix size.factor = 10000 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * 10000),
      as.matrix(GTestimate::GTestimate(test_delayed, size.factor = 10000, log1p.transform = F, rescale = T))
    )
  })
  
  test_that("GTestimate.DelayedMatrix size.factor = 10000 log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))) * 10000),
      as.matrix(GTestimate::GTestimate(test_delayed, size.factor = 10000, log1p.transform = T, rescale = T))
    )
  })
  
  test_that("GTestimate.DelayedMatrix size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(GTestimate::GTestimate(test_delayed, size.factor = test_sf, log1p.transform = T, rescale = T))
    )
  })
}

{
  test_that("GTestimate.Seurat size.factor = 1 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      as.matrix(GetAssayData(GTestimate::GTestimate(pbmc_small, size.factor = 1, log1p.transform = F, rescale = T), slot = 'data', assay = 'GTestimate'))
    )
  })
  
  test_that("GTestimate.Seurat size.factor = 10000 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      (t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))) * 10000),
      as.matrix(GetAssayData(GTestimate::GTestimate(pbmc_small, size.factor = 10000, log1p.transform = F, rescale = T), slot = 'data', assay = 'GTestimate'))
    )
  })
  
  test_that("GTestimate.Seurat size.factor = 10000 log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))) * 10000),
      as.matrix(GetAssayData(GTestimate::GTestimate(pbmc_small, size.factor = 10000, log1p.transform = T, rescale = T), slot = 'data', assay = 'GTestimate'))
    )
  })
  
  test_that("GTestimate.Seurat size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(GetAssayData(GTestimate::GTestimate(pbmc_small, size.factor = test_sf, log1p.transform = T, rescale = T), slot = 'data', assay = 'GTestimate'))
    )
  })
  
  test_that("GTestimate.Seurat size.factor = 1 log1p.transform = FALSE rescale = TRUE missing mass calculation", {
    testthat::expect_equal(
      1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(pbmc_small, size.factor = 1, log1p.transform = F, rescale = T)$missing_mass
    )
  })
  
  test_that("GTestimate.Seurat size.factor = 10000 log1p.transform = TRUE rescale = TRUE missing mass calculation", {
    testthat::expect_equal(
      1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(pbmc_small, size.factor = 10000, log1p.transform = T, rescale = T)$missing_mass
    )
  })
  
  test_that("GTestimate.Seurat size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE missing mass calculation", {
    testthat::expect_equal(
      1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(pbmc_small, size.factor = test_sf, log1p.transform = T, rescale = T)$missing_mass
    )
  })
  
  test_that("GTestimate.Seurat dense size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE missing mass calculation", {
    testthat::expect_equal(
      1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)),
      GTestimate::GTestimate(test_seurat_dense, size.factor = test_sf, log1p.transform = T, rescale = T)$missing_mass
    )
  })
}
{
  test_that("GTestimate.SingleCellExperiment size.factor = 1 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      as.matrix(GTestimate::GTestimate(test_sce, size.factor = 1, log1p.transform = F, rescale = T)@assays@data$GTestimate)
    )
  })
  
  test_that("GTestimate.SingleCellExperiment size.factor = 10000 log1p.transform = FALSE rescale = TRUE", {
    testthat::expect_equal(
      t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))) * 10000 ,
      as.matrix(assay(GTestimate::GTestimate(test_sce, size.factor = 10000, log1p.transform = F, rescale = T), 'GTestimate'))
    )
  })
  
  test_that("GTestimate.SingleCellExperiment size.factor = 10000 log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p((t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))) * 10000)),
      as.matrix(assay(GTestimate::GTestimate(test_sce, size.factor = 10000, log1p.transform = T, rescale = T), 'GTestimate'))
    )
  })
  
  test_that("GTestimate.SingleCellExperiment size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(assay(GTestimate::GTestimate(test_sce, size.factor = test_sf, log1p.transform = T, rescale = T), 'GTestimate'))
    )
  })
  
  test_that("GTestimate.SingleCellExperiment delayed size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE", {
    testthat::expect_equal(
      log1p(t(t(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))/colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0)) * colSums(test_matrix)/test_sf)),
      as.matrix(assay(GTestimate::GTestimate(test_sce_delayed, size.factor = test_sf, log1p.transform = T, rescale = T), 'GTestimate'))
    )
  })
  
  test_that("GTestimate.SingleCellExperiment size.factor = 1 log1p.transform = FALSE rescale = TRUE missing mass calculation", {
    testthat::expect_equal(
      as.numeric(1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      GTestimate::GTestimate(test_sce, size.factor = 1, log1p.transform = F, rescale = T)$missing_mass
    )
  })
  
  test_that("GTestimate.SingleCellExperiment size.factor = 10000 log1p.transform = TRUE rescale = TRUE missing mass calculation", {
    testthat::expect_equal(
      as.numeric(1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      GTestimate::GTestimate(test_sce, size.factor = 10000, log1p.transform = T, rescale = T)$missing_mass
    )
  })
  
  test_that("GTestimate.SingleCellExperiment size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE missing mass calculation", {
    testthat::expect_equal(
      as.numeric(1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      GTestimate::GTestimate(test_sce, size.factor = test_sf, log1p.transform = T, rescale = T)$missing_mass
    )
  })
  
  test_that("GTestimate.SingleCellExperiment delayed size.factor = computePooledFactors() log1p.transform = TRUE rescale = TRUE missing mass calculation", {
    testthat::expect_equal(
      as.numeric(1 - colSums(edgeR::goodTuringProportions(test_matrix) %>% replace(test_matrix==0,0))),
      GTestimate::GTestimate(test_sce_delayed, size.factor = test_sf, log1p.transform = T, rescale = T)$missing_mass
    )
  })
}

test_that("Error messages", {
  testthat::expect_error(
    GTestimate::GTestimate(test_matrix, size.factor = 'A'),
    'Argument size.factor has to be numeric'
  )

  testthat::expect_error(
    GTestimate::GTestimate(test_matrix, log1p.transform = 'A'),
    'Argument log1p.transform has to be logical'
  )

  testthat::expect_error(
    GTestimate::GTestimate(test_matrix, size.factor = 1:10),
    'Length of size.factor must be 1 or equal to number of columns in object'
  )

  testthat::expect_error(
    GTestimate::GTestimate(test_matrix+0.1),
    'Count matrix may only contain integers'
  )
})
