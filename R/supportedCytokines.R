#' Gene set to score for the CytoSig or the Reactome databases.
#' 
#' Returns the names of the cytokines supported by either the CytoSig or the Reactome databases. 
#' 
#' @param database Database used for gene set construction and set scoring.
#' * "cytosig" returns the 41 cytokines scored using the CytoSig database. 
#' * "reactome" returns the 30 cytokines scored using the Reactome database. 
#' @return List of cytokines associated with the CytoSig or the Reactome databases. 
#' @examples
#' supportedCytokines(database = "cytosig")
#' supportedCytokines(database = "reactome")
#'
#' @export
supportedCytokines <- function(database = "cytosig") { 
  if (database == "cytosig") {
    cytokine.list <- c("IL13", "IL4", "NO", "TGFB3", "GDF11", "TGFB1", "MCSF", "GCSF", 
                       "IL10", "WNT3A", "ActivinA", "BMP4", "BMP2", "BMP6", "IL22", 
                       "IFNG", "IL27", "IFNL", "LIF", "IL6", "OSM", "IL12", "IL21",
                       "LTA", "IL3", "IL2", "GMCSF", "IL15", "TRAIL", "TNFSF12", 
                       "IL1A", "CD40LG", "IL1B", "TNFA", "IL17A", "HGF", "BDNF", 
                       "EGF", "VEGFA", "CXCL12", "FGF2")
    return(cytokine.list[order(cytokine.list)])
  }
    else {
      cytokine.list <- c("IL13", "IL4", "TGFB3", "TGFB1", 
                         "IL10", "WNT3A", "BMP4", "BMP2", "BMP6", "IL22", 
                         "IFNG", "IL27", "LIF", "IL6", "OSM", "IL21",
                         "IL3", "IL2", "IL15", "TNFSF12", 
                         "IL1A", "CD40LG", "IL1B", "IL17A", "HGF", "BDNF", 
                         "EGF", "VEGFA", "CXCL12", "FGF2")
      return(cytokine.list[order(cytokine.list)])
    }
}
