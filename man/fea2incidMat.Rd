\name{fea2incidMat}
\alias{fea2incidMat}
\alias{toMatrix}
\alias{adjMatrix}
\title{
Transforms FEA output into incidence matrices.
}
\description{
Transforms the Funtional Enrichment Analysis (FEA) results into cluster-gene incidence matrices. 
}
\usage{
fea2incidMat(feaResults, key = "Genes", sepChar = NULL, clusterColumn = NULL, filterAttribute = NULL, filterOperator = "<", filterThreshold = 0, removeFilteredGtl = NULL)
}
\arguments{
  \item{feaResults}{
list or data.frame/matrix. Output from one of the FEA functions.
}
  \item{key}{
"Genes" or "Terms". To build gene- or term-based networks.
}
  \item{sepChar}{
character. Character separating genes or terms in the same field. By default: "," for genes and ";" for terms.
}
  \item{clusterColumn}{
character. Name of the column that contains the value to group gene-term sets. Only required if it is different than "Cluster" or "Metagroup".
}
  \item{filterAttribute}{
character or data.frame. Attribute to filter the clusters/metagroups. Filtered clusters/metagroups will not be included in the matrices (and subsequent networks).
Its value should be the data.frame column to use for filtering. It can be provided as character (column name) or data.frame (subset of the data.frame with drop=FALSE).
}
  \item{filterOperator}{
character. Logical operator used for filtering. i.e. ">" (bigger than), "<=" (smaller or equal than), "==" (equal), "!=" (different), "\%\%" (included in),... 
The evaluation order is left to right: filterAttribute ">" filterThreshold, will filter out clusters with filter attribute bigger than the threshold.
}
  \item{filterThreshold}{
numeric. Sets the value to compare to.
}
  \item{removeFilteredGtl}{
logical. Only used by GeneTerm Linker term network. If FALSE, it includes generic terms filtered by GeneTerm Linker from final metagroups. 
}
}
\value{
List: 
\item{clustersMatrix or metagroupsMatrix}{Incidende matrix with the genes or Terms in each cluster or metagroup.}
\item{gtSetsMatrix }{Incidende matrix with the genes or Terms in each gene-term set}
\item{filteredOut }{Clusters or metagroups which where filtered out and therefore not included in the incidence matrices. NULL if none.}

}
\seealso{
Next step in the workflow: 
	\code{\link[FGNet:functionalNetwork]{functionalNetwork()}}
	
%%Previous step in the workflow: FEA

Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\examples{
jobID <- 3907019
results <- fea_gtLinker_getResults(jobID)
incidMat <- fea2incidMat(results)

# Filtering (threshold)
incidMat <- fea2incidMat(results, filterAttribute="Silhouette Width", filterThreshold=0.2)

incidMat$filteredOut
head(incidMat$metagroupsMatrix)
head(incidMat$gtSetsMatrix)

functionalNetwork(incidMat)
    
# Term-based network
incidMatTerms <- fea2incidMat(results, key="Terms")
functionalNetwork(incidMatTerms, plotOutput="dynamic")

# Including generic terms filterd by GtLinker from final metagroups:
incidMatTerms <- fea2incidMat(results, key="Terms",removeFilteredGtl=FALSE)
functionalNetwork(incidMatTerms, plotOutput="dynamic", plotType="bipartite")

# Filtering by keyword
keywords <- c("rna")
selectedGroups <- sapply(getTerms(results),
    function(x) any(grep(paste("(", paste(keywords, collapse="|") ,")",sep=""), tolower(x))))

resultsCbind <- results
resultsCbind$metagroups <- cbind(results$metagroups,
    selectedKeywords=as.numeric(selectedGroups))

matSelectedGroups <- fea2incidMat(resultsCbind,
    filterAttribute="selectedKeywords", filterThreshold=1)

functionalNetwork(matSelectedGroups)
getTerms(results)[selectedGroups]
}
