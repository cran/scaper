#' Cytokine activity scores for a Seurat matrix.
#' 
#' Computes cell-level estimates of cytokine activity for a scRNA-seq Seurat
#' count matrix using the scapeForSeurat method. SCAPE activity estimates are computed by scoring weighted gene 
#' sets from the CytoSig or Reactome databases using the Variance-adjusted Mahalanobis (VAM) method 
#' as implemented in the \code{\link[VAM:vamForSeurat]{VAM::vamForSeurat()}} function. Individual gene sets for subsequent scoring can be reconstructed using 
#' the \code{\link[scaper]{genesetCytoSig}} and the \code{\link[scaper]{genesetReactome}} functions 
#' for the CytoSig and the Reactome database, respectively. 
#' @seealso{\code{\link{genesetCytoSig}}, \code{\link{genesetReactome}}, \code{\link{scape}}}
#' 
#' @param seurat.object Seurat counts matrix. 
#' @param database Database used for gene set construction and set scoring.
#' * "cytosig" (default) performs scoring for up to 41 cytokines using the CytoSig database. 
#' * "reactome" performs scoring for up to 30 cytokines using the Reactome database. 
#' @param cytokine Vector of cytokine names to score for activity. The default value of "all" 
#' will score all 41 cytokines supported by CytoSig or 31 supported by Reactome. Please see 
#' function \code{\link[scaper]{supportedCytokines}} to view all the CytoSig or the Reactome 
#' specific scored cytokines. 
#' @param normalize Boolean indicator for whether normalization should be performed before performing gene set scoring. 
#' @return Seurat object consisting of cell-level cytokine activity scores returned as a separate assay (scape for scoring via the CytoSig database and VAMcdf for scoring via the Reactome database). 
#' @examples
#' library(SeuratObject)
#' CytoSig.score.output.all <- scapeForSeurat(seurat.object = pbmc_small, 
#' database = "cytosig", cytokine = "all", normalize=TRUE)
#' (as.data.frame(CytoSig.score.output.all@assays$scape@data))[1:6,1:3]
#' CytoSig.score.output.specific <- scapeForSeurat(seurat.object = pbmc_small, 
#' database = "cytosig", cytokine = c("IL4", "IL13"), normalize=TRUE)
#' (as.data.frame(CytoSig.score.output.specific@assays$scape@data))[,1:3]
#'
#' @export
scapeForSeurat <- function(seurat.object, database="cytosig", cytokine="all", normalize=TRUE) { 
  if (missing(seurat.object)) {
    stop("scRNA-seq Seurat matrix is missing.")
  }
  if (missing(database)) {
    stop("Specify database for gene set scoring.")
  }
  if (normalize == TRUE) {
    seurat.object <- NormalizeData(seurat.object) 
  }
  if (database == "cytosig") {
    if (length(cytokine) == 1 && cytokine == "all") {
      cytokine.list <- supportedCytokines(database = "cytosig")
    } else {
      cytokine.list <- cytokine
    }
    total.output <- mat.cytosig.cytokine.weights
    for (i in 1:length(total.output)) { 
      label.genes <- paste(unique(total.output[[i]]$cytokineLabel), "Genes", sep = "-")
      out.data <- total.output[[i]]
      positive.weighted.data <- out.data %>% dplyr::filter(weight > 0) %>% dplyr::arrange(desc(weight))
      positive.weighted.genes <- positive.weighted.data[positive.weighted.data$gene %in% rownames(seurat.object),]
      positive.weighted.indices <- which(rownames(seurat.object) %in% positive.weighted.genes$gene)
      label.genes.positive <- paste(unique(total.output[[i]]$cytokineLabel), "GenesPositive", sep = "-")
      assign(label.genes.positive, positive.weighted.indices)
      label.weights.positive <- paste(unique(total.output[[i]]$cytokineLabel), "WeightsPositive", sep = "-")
      assign(label.weights.positive, positive.weighted.genes$weight)
      
      negative.weighted.data <- out.data %>% dplyr::filter(weight < 0) %>% dplyr::arrange(desc(weight))
      negative.weighted.genes <- negative.weighted.data[negative.weighted.data$gene %in% rownames(seurat.object),]
      negative.weighted.indices <- which(rownames(seurat.object) %in% negative.weighted.genes$gene)
      label.genes.negative <- paste(unique(total.output[[i]]$cytokineLabel), "GenesNegative", sep = "-")
      assign(label.genes.negative, negative.weighted.indices)
      label.weights.negative <- paste(unique(total.output[[i]]$cytokineLabel), "WeightsNegative", sep = "-")
      assign(label.weights.negative, abs(negative.weighted.genes$weight))
    }
    genes.pattern.positive <- paste(cytokine.list, "-GenesPositive", sep="")
    list_genes_positive <- do.call("list", mget(genes.pattern.positive))
    names(list_genes_positive) <- genes.pattern.positive
    weights.pattern.positive <- paste(cytokine.list, "-WeightsPositive", sep="")
    list_weights_positive <- do.call("list", mget(weights.pattern.positive))
    names(list_weights_positive) <- weights.pattern.positive
    vam.out.positive <- vamForSeurat(seurat.data = seurat.object, 
                                         gene.weights = list_weights_positive,
                                         gene.set.collection = list_genes_positive)
    
    genes.pattern.negative <- paste(cytokine.list, "-GenesNegative", sep="")
    list_genes_negative <- do.call("list", mget(genes.pattern.negative))
    names(list_genes_negative) <- genes.pattern.negative
    weights.pattern.negative <- paste(cytokine.list, "-WeightsNegative", sep="")
    list_weights_negative <- do.call("list", mget(weights.pattern.negative))
    names(list_weights_negative) <- weights.pattern.negative
    vam.out.negative <- vamForSeurat(seurat.data = seurat.object, 
                                         gene.weights = list_weights_negative,
                                         gene.set.collection = list_genes_negative)
    scape.out <- (vam.out.positive@assays$VAMcdf@data + (1-vam.out.negative@assays$VAMcdf@data))/2
    vam.out.positive[['VAMcdf']] <- NULL
    scape.assay <- CreateAssayObject(scape.out)
    vam.out.positive[['scape']] <- scape.assay
    rownames(vam.out.positive@assays$scape@data) <- cytokine.list
    return(vam.out.positive)
  }
  else {
    if (cytokine == "all") {
      cytokine.list <- supportedCytokines(database = "reactome")
    } else {
      cytokine.list <- cytokine
    }
    total.output <- mat.reactome.cytokine
    for (i in 1:length(total.output)) { 
      label.genes <- paste(unique(total.output[[i]]$cytokineLabel), "Genes", sep = "-")
      out.data <- total.output[[i]]
      out.data.genes <- out.data[out.data$gene %in% rownames(seurat.object),]
      out.data.indices <- which(rownames(seurat.object) %in% out.data.genes$gene)
      label.genes <- paste(unique(total.output[[i]]$cytokineLabel), "Genes", sep = "-")
      assign(label.genes, out.data.indices)
    }
    genes.pattern <- paste(cytokine.list, "-Genes", sep="")
    genes.pattern.list <- do.call("list", mget(genes.pattern))
    list_genes <- genes.pattern.list
    names(list_genes) <- genes.pattern
    list_genes2 <- list_genes[lengths(list_genes) >  0]
    vam.out <- vamForSeurat(seurat.data = seurat.object, 
                            gene.set.collection = list_genes2)
    rownames(vam.out@assays$VAMcdf@data) <- gsub("-Genes", "", 
                                                 rownames(vam.out@assays$VAMcdf@data))
    return(vam.out)
  } 
}