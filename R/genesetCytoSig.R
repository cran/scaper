#' CytoSig gene set construction.
#' 
#' Returns the gene set associated with input cytokine(s) from the CytoSig database
#' given manually specified cytokine-specific output csv file(s) under 
#' the extdata directory with file name beginning with 
#' the specified cytokine (i.e., 'IL6_output.csv') as currently provided for the IL6 cytokine. 
#' 
#' @param cytokine.eval Cytokine(s) associated with the query.
#' @param file.name List of XML file(s) associated with the cytokine(s) beginning with the specific cytokine name.
#' @return Dataframe consisting of genes and the associated log2 fold change values associated with the specific cytokine(s). 
#' @examples
#' file.name.cytosig1 <- system.file("extdata", "IL6_output.csv", package = "scaper")
#' genesetCytoSig(cytokine.eval = "IL6", file.name = file.name.cytosig1) %>% head(10)
#' @export
genesetCytoSig <- function(cytokine.eval, file.name) {
  total.output <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(total.output) <- c("gene", "weight", "cytokine")
  for (j in 1:length(cytokine.eval)) {
    file.name.base <- basename(file.name) 
    path.cytokine <- file.name.base[grepl(paste("^",  cytokine.eval[j],".*", ".csv$", sep=""), file.name.base)]
    if (is.null(path.cytokine) || length(path.cytokine) == 0) {
      stop("Specify cytokine-specific output csv file from the CytoSig database with file name beginning with the specified cytokine.")
    }
    path.out <- read.csv(file.name[grepl(path.cytokine, file.name)])
    path.out$Value <- as.numeric(as.character(path.out$Value))
    path.out2 <- path.out %>%
      dplyr::arrange(desc(Value)) %>%
      dplyr::group_by(Gene) %>%
      dplyr::summarise(meanVal = mean(Value))
    total.genes <- path.out2$Gene
    weights.genes <- path.out2$meanVal
    total.out <- as.data.frame(cbind(total.genes, weights.genes, cytokine.eval[j]))
    total.out <- total.out %>% dplyr::arrange(desc(weights.genes))
    names(total.out) <- c("gene", "weight", "cytokine")
    total.output <- rbind(total.output, total.out)
    total.output <- total.output %>% dplyr::arrange(cytokine)
  }
  return(total.output)
}
