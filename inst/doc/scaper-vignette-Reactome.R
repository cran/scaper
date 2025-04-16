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
genesetReactome(cytokine.eval = "IL6", 
                file.name = system.file("extdata", 
                                        "IL6_Interleukin6_signaling.xml", 
                                        package = "scaper")) %>% head(10)
supportedCytokines(database = "reactome")

## ----message=FALSE, warning=FALSE---------------------------------------------
Reactome.score.output <- scapeForSeurat(seurat.object = pbmc_small, 
                                       database = "reactome", 
                                       normalize = TRUE)
reactome_mat <- as.data.frame(t(as.matrix(Reactome.score.output@assays$VAMcdf@data)))
pheatmap(reactome_mat, fontsize_row = 4, fontsize_col = 7, 
         cluster_rows = FALSE, cluster_cols = FALSE)

## ----message=FALSE, warning=FALSE---------------------------------------------
pheatmap(reactome_mat, fontsize_row = 4, fontsize_col = 7)

