\name{toMatrix}
\alias{toMatrix}
\alias{adjMatrix}
\title{
Transforms into group-genes incidence matrices.
}
\description{
Transforms the raw results from GeneTermLinker or DAVID into incidence matrices. 
}
\usage{
toMatrix(results, attribute=NULL, threshold=0, 
key="Genes", removeFiltered=NULL)
}
\arguments{
  \item{results}{
list or data.frame. \emph{results} returned by \code{\link[FGNet:getResults_gtLinker]{getResults_gtLinker}} or \code{\link[FGNet:getResults_david]{getResults_david}}.
}
  \item{attribute}{
data.frame. Attribute (column from results) to filter the metagroups/clusters.
}
  \item{threshold}{
numeric. Metagroups/cluster with \emph{attribute} lower than this threshold will be filtered.
}
  \item{key}{
"Genes" or "Terms". To build genes- or terms-based networks.
}
  \item{removeFiltered}{
logical. If FALSE, it includes generic terms filtered by GeneTerm Linker from final metagroups. Only available when building a terms network with GeneTerm Linker.
}
}
\value{
List: 
\item{metagroupsMatrix or clustersMatrix }{Incidende matrix: Genes or Terms in each metagroup or cluster.}
\item{gtSetsMatrix }{Incidende matrix: Genes or Terms in each gene-term set}
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
incidMat <- toMatrix(results)

# Filtering (threshold)
incidMat <- toMatrix(results, 
 attribute=results$metagroups[,"Silhouette Width", drop=FALSE], threshold=0.2)

incidMat$filteredOut
head(incidMat$metagroupsMatrix)
head(incidMat$gtSetsMatrix)

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

# Term-based network
incidMatTerms <- toMatrix(results, key="Terms")
functionalNetwork(incidMatTerms, plotType="dynamic")

# Including generic terms filterd by GtLinker from final metagroups:
incidMatTerms2 <- toMatrix(results, key="Terms", removeFiltered=FALSE)
intersectionNetwork(incidMatTerms2)
}
