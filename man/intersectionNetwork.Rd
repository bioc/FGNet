\name{intersectionNetwork}
\alias{intersectionNetwork}
\title{
Metagroup/Cluster intersection network.
}
\description{
Plots a simplified version of the functional network, only with the genes in more than one metagroup/cluster.
}
\usage{
intersectionNetwork(metagroupsMatrix, plotType = "dynamic", 
vLayout = "kk", returnGraph = FALSE, vSize = 12, vLabelCex = 2/3, 
legendMg=NULL, grPrefix="", plotTitle = "Nodes in several metagroups",
keepColors=TRUE, plotAllMg=FALSE)
}
\arguments{
  \item{metagroupsMatrix}{
list or matrix. raw output (list) from \code{\link[FGNet:toMatrix]{toMatrix}}, or 'metagroupsMatrix' or 'clustersMatrix'.
}
  \item{plotType}{
character "static", "dynamic" or "none". "static" will generate a standard plot. "dynamic" will produce an interactive tkplot (metagroups background cannot be drawn). "none" will not plot the network.
}
  \item{vLayout}{
character.  "kk" (Kamada Kawai), "circle", or "sugiyama" (hierarquical).
}
  \item{returnGraph}{
logical. If TRUE, the igraph object containing the network is returned.
}
  \item{vSize}{
numeric. Vertex size.
}
  \item{vLabelCex}{
numeric. A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default label size. 
}
  \item{legendMg}{
character. Label to show next to the metagroup/cluster id in the node label.
}
  \item{grPrefix}{
character. Prefix for the metagroup/cluster.
}
  \item{plotTitle}{
character. Title to show on the plot.
}
  \item{keepColors}{
logical. If TRUE, it will keep the same colors for all the plots, independently of the filtered groups. Only available if metagroupGenesMatrix is the raw result from \code{\link[FGNet:toMatrix]{toMatrix}}.
}
  \item{plotAllMg}{
logical. If TRUE, metagroup/cluster nodes will be plotted even if they are not connected to any more nodes.
}
}
\value{
 Plots the network.
 If returnGraph=TRUE: Returns the igraph object with the network.
}

\seealso{
	Full description of the package:
		\code{\link[FGNet:FGNet-package]{FGNet}}
}
\examples{
jobID <- 3907019
results <- getResults_gtLinker(jobID, jobName="gtLinkerExample")
incidMat <- toMatrix(results, 
 attribute=results$metagroups[,"Silhouette Width", drop=FALSE], threshold=0)

intersectionNetwork(incidMat)
intNw <- intersectionNetwork(incidMat, 
 vLayout="sugiyama", plotType="static", returnGraph=TRUE) 
 
 
# Term-based network
incidMatTerms <- toMatrix(results, key="Terms")
intersectionNetwork(incidMatTerms)

# Including generic terms filterd by GtLinker from final metagroups:
incidMatTerms2 <- toMatrix(results, key="Terms", removeFiltered=FALSE)
intersectionNetwork(incidMatTerms2)
}
