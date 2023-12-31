% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genesetReactome.R
\name{genesetReactome}
\alias{genesetReactome}
\title{Reactome gene set construction.}
\usage{
genesetReactome(cytokine.eval, file.name)
}
\arguments{
\item{cytokine.eval}{Cytokine(s) associated with the query.}

\item{file.name}{List of XML file(s) associated with the cytokine(s) beginning with the specific cytokine name.}
}
\value{
Dataframe consisting of list of genes on the molecular pathway associated with the specified cytokine(s).
}
\description{
Returns the gene set consisting of genes on the Reactome pathway
hierarchy given input cytokine(s) and specified list of pathway
hierarchy xml file(s) with file name beginning with the specified cytokine
(i.e., 'IL6_Interleukin6_signaling.xml') as currently provided for the IL6 cytokine.
}
\examples{
file.name.reactome1 <- system.file("extdata", "IL6_Interleukin6_signaling.xml", 
package = "scaper")
genesetReactome(cytokine.eval = "IL6", file.name = file.name.reactome1) \%>\% head(10)
}
