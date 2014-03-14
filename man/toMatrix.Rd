\name{toMatrix}
\alias{toMatrix}
\alias{adjMatrix}
\title{
Transforms into group-genes incidence matrices.
}
\description{
Transforms the raw results from GeneTermLinker or DAVID into genes-metagroups and genes-gtsets incidence matrices. 
}
\usage{
toMatrix(geneTermSets, attribute = NULL, threshold = 0)
}
\arguments{
  \item{geneTermSets}{
data.frame. \emph{geneTermSets} returned by \code{\link[FGNet:getResults_gtLinker]{getResults_gtLinker}} or \emph{clusters} returned by \code{\link[FGNet:getResults_david]{getResults_david}}.
}
  \item{attribute}{
data.frame. Attribute (column from results) to filter the metagroups/clusters.
}
  \item{threshold}{
numeric. Metagroups/cluster with \emph{attribute} lower than this threshold will be filtered.
}
}
\value{
List: 
\item{metagroupGenesMatrix }{Table 'Genes - Metagroups' or 'Genes - Clusters'.}
\item{gtSetGenesMatrix }{Table 'Genes - Gene-Term sets'.}
\item{filteredOut }{Metagroups/clusters which where filtered out and therefore not included in the incidence matrices. NULL if none.}
}
\seealso{
	Next step in the workflow: 
		\code{\link[FGNet:functionalNetwork]{functionalNetwork()}} (Plots)\cr
		
	Previous step in the workflow: 
		\code{\link[FGNet:getResults_gtLinker]{getResults_gtLinker()}} or \code{\link[FGNet:getResults_david]{getResults_david()}}\cr
		
	Full description of the package:
		\code{\link[FGNet:FGNet-package]{FGNet}}
}
\examples{
jobID <- 3907019
results <- getResults_gtLinker(jobID)
incidMat <- toMatrix(results$geneTermSets)

# Filtering (threshold)
incidMat <- toMatrix(results$geneTermSets, 
 attribute=results$metagroups[,"Silhouette Width", drop=FALSE], threshold=0.2)

incidMat$filteredOut
head(incidMat$metagroupGenesMatrix)
head(incidMat$gtSetGenesMatrix)

functionalNetwork(incidMat)

# Filtering (keyword)
keywords <- c("rna")
selectedGroups <- sapply(getTerms(results),
function(x)
any(grep(paste("(", paste(keywords, collapse="|") ,")",sep=""), tolower(x))))

resultsCbind <- results
resultsCbind$metagroups <- cbind(results$metagroups,
 selectedKeywords=as.numeric(selectedGroups))
 
matSelectedGroups <- toMatrix(resultsCbind$geneTermSets,
 attribute=resultsCbind$metagroups[,"selectedKeywords", drop=FALSE], threshold=1)
 
functionalNetwork(matSelectedGroups)
}
