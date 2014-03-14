\name{query_david}
\alias{query_david}
\title{
Query DAVID.
}
\description{
Sends a new analysis to DAVID [1] to perform the functional enrichment and clustering.
}
\usage{
query_david(genes, geneIdType = "ENSEMBL_GENE_ID", 
 annotations = c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", 
 "KEGG_PATHWAY", "INTERPRO"), email=NULL, 
 argsWS = c(overlap=4L, initialSeed=4L, finalSeed=4L, linkage=0.5, kappa=35L))
}
\arguments{
  \item{genes}{
character vector. List of genes to analyze.
}
  \item{geneIdType}{
character vector. Type of gene identifier.\cr
Web API: ENSEMBL_GENE_ID, ENTREZ_GENE_ID, GENE_SYMBOL, UNIPROT_ID... For more, check DAVID's API documentation.\cr
Web Service: run \code{getIdTypes(DAVIDWebService$new(email=...))}
}
  \item{annotations}{
character vector. Annotation spaces for the functional analysis. \cr
Web API: check DAVID's API documentation.\cr
Web Service: run \code{getAllAnnotationCategoryNames(DAVIDWebService$new(email=...))}.
}
  \item{email}{
character. If provided, the query will be performed though DAVID's Web Service (recommended). Requires registration.
}
  \item{argsWS}{
named integer vector. Additional arguments for the clustering. Only available using the web service. 
}
}
\value{
Returns the url/location of the text file containing the results of the analysis.
}
\details{

If an email is provided, the query will be performed through the web service, if not through the API:


\bold{Web service:} Recommended option. Requiers registration at \url{http://david.abcc.ncifcrf.gov/webservice/register.htm}.

\bold{API:} Allows to perform a small query without registering. \cr
Maximum number of genes: 400. It can be less depending on the ID types.\cr
The default arguments for the clustering performed through the API are:
\code{web_API_defaults <- c(overlap=3L, initialSeed=3L, finalSeed=3L, linkage=0.5, kappa=50L)}\cr
More details and full list of gene ID types and annotations are available at: 
\url{http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html}.\cr

\bold{Website:} If the \emph{functional annotation and clustering} has been performed directly at DAVID's website (\url{http://david.abcc.ncifcrf.gov/summary.jsp}) this function is not required. Provide the file (or the URL of the file) containing the results of the analysis to \code{\link[FGNet:getResults_david]{getResults_david()}}.
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
\dontrun{
genesYeast <- c("YBL084C", "YDL008W", "YDR118W", "YDR301W", "YDR448W", 
"YFR036W", "YGL240W", "YHR166C", "YKL022C", "YLR102C", "YLR115W", "YLR127C",
"YNL172W", "YOL149W", "YOR249C")

# Default:
txtFile <- query_david(genesYeast)

# Web service (requires email):
txtFile_ws <- query_david(genesYeast, email="example@email.com")

# The API has different default arguments than the web service.
# To obtain exactly the same results using the web service, use:
API_defaults <- c(overlap=3L, initialSeed=3L, finalSeed=3L, linkage=0.5, kappa=50L) 
query_david(genesYeast, email="...", argsWS=API_defaults)
}
}
