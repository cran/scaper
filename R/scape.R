#' Cytokine activity scores for a normalized matrix.
#' 
#' Computes cell-level estimates of cytokine activity for a normalized scRNA-seq 
#' count matrix using the SCAPE method. SCAPE activity estimates are computed by scoring weighted genes 
#' sets from the CytoSig or Reactome databases using the \code{\link[VAM:vamForCollection]{VAM::vamForCollection()}}
#' function. Individual gene sets for subsequent scoring can be reconstructed using 
#' the \code{\link[scaper]{genesetCytoSig}} and the \code{\link[scaper]{genesetReactome}} functions 
#' for the CytoSig and the Reactome database, respectively. 
#' @seealso{\code{\link{genesetCytoSig}}, \code{\link{genesetReactome}}}
#' 
#' @param counts.matrix A \eqn{m x n} normalized counts matrix with \eqn{m} samples and \eqn{n} genes. 
#' @param database Database used for gene set construction and set scoring.
#' * "cytosig" performs scoring for up to 41 cytokines using the CytoSig database. 
#' * "reactome" performs scoring for up to 30 cytokines using the Reactome database. 
#' @param cytokine Vector of cytokine names to score for activity. The default value 
#' of "all" will score all 41 cytokines supported by CytoSig or 31 supported by Reactome. Please see 
#' function \code{\link[scaper]{supportedCytokines}} to view all the CytoSig or the Reactome 
#' specific scored cytokines. 
#' @return A \eqn{m x p} matrix consisting of the cell-level cytokine activity scores for \eqn{p} cytokines. 
#' @examples
#' library(Seurat)
#' library(SeuratObject)
#' pbmc_small <- NormalizeData(pbmc_small)
#' counts.matrix <- as.data.frame(t(as.matrix(pbmc_small@assays$RNA@data)))
#' CytoSig.score.output <- scape(counts.matrix = counts.matrix, 
#' database = "cytosig")
#' head(CytoSig.score.output)[,1:3]
#' CytoSig.score.output.specific <- scape(counts.matrix = counts.matrix, 
#' database = "cytosig", cytokine = c("IL4", "IL13"))
#' head(CytoSig.score.output.specific)
#'
#' @export
scape <- function(counts.matrix, database="cytosig", cytokine="all") { 
  if (missing(counts.matrix)) {
    stop("scRNA-seq counts matrix is missing.")
  }
  if (database == "cytosig") {
    if (length(cytokine) == 1 && cytokine == "all") {
      cytokine.list <- supportedCytokines(database = "cytosig")
    } else {
      cytokine.list <- cytokine
    }
    total.output <- mat.cytosig.cytokine.weights[cytokine.list]
    for (i in 1:length(total.output)) {
      label.genes <- paste(unique(total.output[[i]]$cytokineLabel), "Genes", sep = "-")
      out.data <- total.output[[i]]
      positive.weighted.data <- out.data %>% dplyr::filter(weight > 0) %>% dplyr::arrange(desc(weight))
      positive.weighted.genes <- positive.weighted.data[positive.weighted.data$gene %in% colnames(counts.matrix),]
      positive.weighted.indices <- which(colnames(counts.matrix) %in% positive.weighted.genes$gene)
      label.genes.positive <- paste(unique(total.output[[i]]$cytokineLabel), "GenesPositive", sep = "-")
      assign(label.genes.positive, positive.weighted.indices)
      label.weights.positive <- paste(unique(total.output[[i]]$cytokineLabel), "WeightsPositive", sep = "-")
      assign(label.weights.positive, positive.weighted.genes$weight)
      
      negative.weighted.data <- out.data %>% dplyr::filter(weight < 0) %>% dplyr::arrange(desc(weight))
      negative.weighted.genes <- negative.weighted.data[negative.weighted.data$gene %in% colnames(counts.matrix),]
      negative.weighted.indices <- which(colnames(counts.matrix) %in% negative.weighted.genes$gene)
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
    vam.out.positive <- vamForCollection(gene.expr = counts.matrix, 
                            gene.weights = list_weights_positive,
                            gene.set.collection = list_genes_positive)
    
    genes.pattern.negative <- paste(cytokine.list, "-GenesNegative", sep="")
    list_genes_negative <- do.call("list", mget(genes.pattern.negative))
    names(list_genes_negative) <- genes.pattern.negative
    weights.pattern.negative <- paste(cytokine.list, "-WeightsNegative", sep="")
    list_weights_negative <- do.call("list", mget(weights.pattern.negative))
    names(list_weights_negative) <- weights.pattern.negative
    vam.out.negative <- vamForCollection(gene.expr = counts.matrix, 
                                         gene.weights = list_weights_negative,
                                         gene.set.collection = list_genes_negative)
    vam.out.positive$cdf.value <- (vam.out.positive$cdf.value + (1-vam.out.negative$cdf.value))/2
    colnames(vam.out.positive$cdf.value) <- cytokine.list
    return(vam.out.positive$cdf.value)
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
      out.data.genes <- out.data[out.data$gene %in% colnames(counts.matrix),]
      out.data.indices <- which(colnames(counts.matrix) %in% out.data.genes$gene)
      label.genes <- paste(unique(total.output[[i]]$cytokineLabel), "Genes", sep = "-")
      assign(label.genes, out.data.indices)
    }
    genes.pattern <- paste(cytokine.list, "-Genes", sep="")
    genes.pattern.list <- do.call("list", mget(genes.pattern))
    list_genes <- genes.pattern.list
    names(list_genes) <- genes.pattern
    list_genes2 <- list_genes[lengths(list_genes) >  0]
    vam.out <- vamForCollection(gene.expr = counts.matrix, 
                            gene.set.collection = list_genes2)
    colnames(vam.out$cdf.value) <- gsub("-Genes", "", 
                                                 colnames(vam.out$cdf.value))
    return(vam.out$cdf.value)
  } 
}