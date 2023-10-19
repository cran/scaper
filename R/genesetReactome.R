#' Reactome gene set construction.
#' 
#' Returns the gene set consisting of genes on the Reactome pathway 
#' hierarchy given input cytokine(s) and specified list of pathway 
#' hierarchy xml file(s) with file name beginning with the specified cytokine 
#' (i.e., 'IL6_Interleukin6_signaling.xml') as currently provided for the IL6 cytokine. 
#' 
#' @param cytokine.eval Cytokine(s) associated with the query.
#' @param file.name List of XML file(s) associated with the cytokine(s) beginning with the specific cytokine name. 
#' @return Dataframe consisting of list of genes on the molecular pathway associated with the specified cytokine(s). 
#' @examples
#' file.name.reactome1 <- system.file("extdata", "IL6_Interleukin6_signaling.xml", 
#' package = "scaper")
#' genesetReactome(cytokine.eval = "IL6", file.name = file.name.reactome1) %>% head(10)
#' @export
genesetReactome <- function(cytokine.eval, file.name) {
  str.output.final.parsed.total <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(str.output.final.parsed.total) <- c("gene", "cytokine")
  for (j in 1:length(cytokine.eval)){
    file.name.base <- basename(file.name) 
    filenames <- file.name.base[grepl(paste("^",  cytokine.eval[j],".*", ".xml$", sep=""), file.name.base)]
    if (is.null(filenames) || length(filenames) == 0) {
      stop("Specify Reactome pathway hierarchy xml file(s) with file name beginning with the specified cytokine.")
    }
    filenames.update <- c()
    for (k in 1:length(filenames)) {
      filenames.update[k] <- file.name[grepl(filenames[k], file.name)]
    }
    path.out <- lapply(filenames.update, read_xml)
    str.output.final <- list()
    for (i in 1:length(path.out)) {
      allPathwayComponents <- xml_find_all(path.out[[i]], "/rdf:RDF/bp:biochemicalReaction") %>% xml_path()
      downstreamSteps <- c()
      for (r in 1:length(allPathwayComponents)) {
        stringInquiry <- paste("/rdf:RDF/bp:biochemicalReaction[", r, "]/bp:NAME", sep = "")
        downstreamSteps[r] <- xml_text(xml_find_all(path.out[[i]], stringInquiry))
      }
      downstreamSteps2 <- unlist(downstreamSteps)
      out <- strsplit(downstreamSteps2, "\\s+") 
      out2 <- strsplit(unlist(out), ":")
      str.output <- str_trim(str_extract(unlist(out2), "([[:upper:][0-9]]){3,}"))
      str.output2 <- unique(str.output[!is.na(str.output)])
      str.output.final[[i]] <- str.output2[grepl("\\D", str.output2)]
    }
    str.output.final.parsed <- as.data.frame(cbind(unlist(str.output.final), 
                                                   cytokine.eval[j]))
    names(str.output.final.parsed) <- c("gene", "cytokine")
    receptor.sub <- ramilowski.database %>% 
      dplyr::filter(Ligand.ApprovedSymbol == cytokine.eval[j])
    receptor.sub1 <- receptor.sub %>% dplyr::select(Receptor.ApprovedSymbol, Ligand.ApprovedSymbol)
    names(receptor.sub1) <- c("gene", "cytokine")
    str.output.final.parsed <- rbind(str.output.final.parsed, receptor.sub1)
    str.output.final.parsed.total <- unique(rbind(str.output.final.parsed.total, str.output.final.parsed))
    str.output.final.parsed.total <- str.output.final.parsed.total %>% dplyr::arrange(cytokine)
  }
  return(str.output.final.parsed.total)
}
