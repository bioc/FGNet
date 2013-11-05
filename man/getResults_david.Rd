\name{getResults_david}
\alias{getResults_david}
\title{
Get results from an analysis with David.
}
\description{
Retrieves the clustering results from David. The .txt file URL can be obtained from \link[FGNet:query_david]{query_david} or from a analysis performed directly at David (http://david.abcc.ncifcrf.gov/summary.jsp) . 
}
\usage{
getResults_david(inputFileLocation, path = getwd(), jobName = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputFileLocation}{
	character. URL of the file to download, or local file if already downloaded.
}
  \item{path}{
  character. Directory in which to save the files.
}
  \item{jobName}{
  character. Folder name and prefix for the files.
}
}
\details{
\emph{\link[FGNet:query_david]{query_david()}} uses DAVID's API to perform the query, therefore the maximum number of genes is limited to 400 or less depending on the used ID and final URL length. In order to perform analyses with more genes, the web interface can be used (\url{http://david.abcc.ncifcrf.gov/summary.jsp}). Once the 'functional annotation and clustering' is done though the website, to continue the workflow through R just provide as 'inputFileLocation' the URL of the .txt downloadable file.
}
\value{
List:\cr
\code{clusters } data.frame containing the clusters and their information:
   \itemize{ 
  \item Cluster: Cluster ID.
  \item EnrichmentScore: Score for the cluster.
  \item nGenes: Number of genes in the cluster.
  \item Genes: Genes in the cluster.
 \item Terms: Terms in the cluster.
 }
\code{geneTermSets } data.frame containing the gene-term sets that support each cluster.
   \itemize{ 
 \item Cluster: Number (id) of the cluster the gene-term set belongs to.
 \item Term: Term in the gene-term set.
 \item Category: Type of annotation of the term (i.e. GO, Kegg...)
 \item Genes: Genes in the gene-term set.
 \item Other stats: Count, PValue, List.Total, Pop.Hits, Pop.Total, Fold.Enrichment, Bonferroni, Benjamini, FDR.
 }
\code{fileName} and location of the raw .txt downloaded from David.
 
}

\seealso{
	Next step in the workflow: 
		\code{\link[FGNet:adjMatrix]{adjMatrix()}}\cr
	
	Previous step in the workflow: 
		\code{\link[FGNet:query_david]{query_david()}}\cr
		
	Equivalent function for GeneTerm Linker:
	\code{\link[FGNet:getResults_gtLinker]{getResults_gtLinker}}\cr
	
	Full description of the package:
	\code{\link[FGNet:FGNet-package]{FGNet}}
}

\examples{
genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9",
 "CDC16", "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1")

txtFile <- query_david(genesYeast)

results <- getResults_david(txtFile)

# If the analysis results are already in a local txt file:
# getResults_david(".../DavidClustering.txt")
}
