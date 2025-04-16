## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(scaper)
library(Seurat)
library(SeuratObject)
library(pheatmap)

## ----message=FALSE, warning=FALSE---------------------------------------------
supportedCytokines(database = "cytosig")
genesetCytoSig(cytokine.eval = "IL6", 
               file.name = system.file("extdata", "IL6_output.csv",
                                       package = "scaper")) %>% head(10)

## ----message=FALSE, warning=FALSE---------------------------------------------
CytoSig.score.output <- scapeForSeurat(seurat.object = pbmc_small, 
                                       database = "cytosig", cytokine = "all", normalize = TRUE)
class(CytoSig.score.output)
GetAssay(object = CytoSig.score.output, assay = "scape")
#DefaultAssay(CytoSig.score.output) <- "scape"
#GetAssayData(object = CytoSig.score.output, slot = "data")
cytosig_mat <- as.data.frame(t(as.matrix(CytoSig.score.output@assays$scape@data)))
pheatmap(cytosig_mat, fontsize_row = 4, fontsize_col = 7, 
         cluster_rows = FALSE, cluster_cols = FALSE)

## ----message=FALSE, warning=FALSE---------------------------------------------
pheatmap(cytosig_mat, fontsize_row = 4, fontsize_col = 7)

## ----message=FALSE, warning=FALSE---------------------------------------------
pbmc_small <- NormalizeData(pbmc_small)
counts.matrix2 <- as.data.frame(t(as.matrix(pbmc_small@assays$RNA@data)))
CytoSig.score.output <- scape(counts.matrix = counts.matrix2, 
                              database = "cytosig", cytokine = "all")
pheatmap(CytoSig.score.output, fontsize_row = 4, fontsize_col = 7, 
         cluster_rows = FALSE, cluster_cols = FALSE)

