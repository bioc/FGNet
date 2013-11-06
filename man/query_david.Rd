\name{query_david}
\alias{query_david}
\title{
Query David.
}
\description{
Sends a new analysis to GeneTerm Linker [1] to perform the functional enrichment and clustering.
}
\usage{
query_david(genes, geneIdType = "GENE_SYMBOL", 
 annotations = c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", 
 "KEGG_PATHWAY", "INTERPRO"))
}
\arguments{
  \item{genes}{
character vector. List of genes to analyze.
}
  \item{geneIdType}{ 
character vector. Type of gene id for the genes provided.
Sample ID types: GENE_SYMBOL, ENSEMBL_GENE_ID, ENTREZ_GENE_ID, UNIPROT_ID. For more, check David's API documentation.
}
  \item{annotations}{
character vector. Annotation spaces for the functional analysis. 
Sample spaces: "GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "INTERPRO". For more, check David's API documentation.
}
}
\value{
Returns the url for the text file containing the results of the analysis. This url which will be required for the next steps.
}
\details{
This function uses David's API. More details and full list of available parameters  (gene ID types and annotations) are available at: 
\url{http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html}. 

This API limits the maximum number of genes to 400 (or less depending on the IDs and the final URL lenght). In order to perform analyses with more genes, just perform the analysis directly at DAVID website (\url{http://david.abcc.ncifcrf.gov/summary.jsp}). Once the 'functional annotation and clustering' is done though the website, to continue the workflow through R just provide as 'inputFileLocation' the URL of the .txt downloadable file.
}
\references{
	[1] Huang DW, Sherman BT, Lempicki RA (2009) Bioinformatics enrichment tools: paths toward the comprehensive functional analysis of large gene lists. Nucleic Acids Res. 37(1):1-13.
}
\seealso{
	Next step in the workflow: 
		\code{\link[FGNet:getResults_david]{getResults_david()}}\cr
		
	Equivalent function for GeneTerm Linker: 
		\code{\link[FGNet:query_gtLinker]{query_gtLinker()}}\cr
		
	Full description of the package:
		\code{\link[FGNet:FGNet-package]{FGNet}}
}
\examples{
genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9",
 "CDC16", "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1")

txtFile <- query_david(genesYeast)
}
