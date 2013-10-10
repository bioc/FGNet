\name{intersectionNetwork}
\alias{intersectionNetwork}
\title{
Metagroup/Cluster intersection network.
}
\description{
Plots a simplified version of the functional network, only with the genes in more than one metagroup/cluster.
}
\usage{
intersectionNetwork(metagroupGenesMatrix, plotType = "dynamic", 
vLayout = "kk", returnGraph = FALSE, vSize = 12, vLabelCex = 2/3, 
legendMg=NULL, grPrefix="", plotTitle = "Nodes in several metagroups")
}
\arguments{
  \item{metagroupGenesMatrix}{
matrix. \emph{metagroupGenesMatrix} returned by  \code{\link[FGNet:adjMatrix]{adjMatrix}}.
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
adjMat <- adjMatrix(results$geneTermSets, 
 attribute=results$metagroups[,"Silhouette Width", drop=FALSE], threshold=0)

intersectionNetwork(adjMat$metagroupGenesMatrix)
intNw <- intersectionNetwork(adjMat$metagroupGenesMatrix, 
 vLayout="sugiyama", plotType="static", returnGraph=TRUE) 
}

