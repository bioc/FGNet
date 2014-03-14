\name{functionalNetwork}
\alias{functionalNetwork}
\title{
Create and plot the functional gene network.
}
\description{
Plots the gene network. Edges join nodes with  common gene-term sets. Background and node color represent metagroup/clusters. White nodes are in several metagroups/clusters.
}
\usage{

functionalNetwork(metagroupGenesMatrix, gtSetGenesMatrix=NULL, 
plotType="static", returnGraph=FALSE,
plotTitle="Functional Network", vSize=12, vLabelCex=3/4, 
vLayout=NULL, bgTransparency=0.4, legendMg=NULL, keepColors=TRUE)
}




\arguments{
 \item{metagroupGenesMatrix}{
matrix or list. \emph{$metagroupGenesMatrix} (matrix) returned by \code{\link[FGNet:toMatrix]{toMatrix}}. If raw output (list) from \code{\link[FGNet:toMatrix]{toMatrix}}, it will be automatically splitted into $metagroupGenesMatrix and $gtSetGenesMatrix.
}
  \item{gtSetGenesMatrix}{
matrix. \emph{$gtSetGenesMatrix} returned by \code{\link[FGNet:toMatrix]{toMatrix}}.
}
  \item{plotType}{
character "static", "dynamic" or "none". "static" will generate a standard plot. "dynamic" will produce an interactive tkplot. In the interactive plot metagroups background cannot be drawn, instead the intersection network will also be shown. "none" will not plot the network.
}
  \item{returnGraph}{
logical. If TRUE, the igraph object containing the network is returned.
}
  \item{plotTitle}{
character. Title to show on the plot.
}
  \item{vSize}{
numeric. Vertex size.
}
  \item{vLabelCex}{
numeric. A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default label size. 
}
  \item{vLayout}{
2 x n matrix. Where n is the number of nodes in the graph, each column gives the (x, y)-coordinates for the corresponding node.
}
  \item{bgTransparency}{
numeric. Value between 0 and 1 for the transparency of the metagroups background.
}
  \item{legendMg}{
character. Label to show next to the metagroup/cluster id in the legend. If FALSE, legend is not shown.
}
  \item{keepColors}{
logical. If TRUE, it will keep the same colors for all the plots, independently of the filtered groups. Only available if metagroupGenesMatrix is the raw result from \code{\link[FGNet:toMatrix]{toMatrix}}.
}
}

\value{
 Plots the functional network.\cr
 If plotType="dynamic" it also plots the \link[FGNet:intersectionNetwork]{intersectionNetwork}.\cr
 If returnGraph=TRUE: Returns the igraph object with the network.\cr
}

\seealso{
  To see the terms included in each metagroup or cluster: 
    \code{\link[FGNet:getTerms]{getTerms()}}
	Full description of the package:
		\code{\link[FGNet:FGNet-package]{FGNet}}
}
\examples{
jobID <- 3907019
results <- getResults_gtLinker(jobID, jobName="gtLinkerExample")
incidMat <- toMatrix(results$geneTermSets, 
										attribute=results$metagroups[,"Silhouette Width", drop=FALSE], threshold=0)

functionalNetwork(incidMat) 
functionalNetwork(incidMat, plotType="dynamic") 
getTerms(results)

# To modify the layout and plot as static network (with metagroup background)...
library(igraph)
saveLayout <- tkplot.getcoords(1)   # tkp.id (ID of the tkplot window)
functionalNetwork(incidMat, vLayout=saveLayout)

# Only return the network, without plotting
fNw <- functionalNetwork(incidMat, plotType="none", returnGraph=TRUE) 
class(fNw)
betweenness(fNw)
igraph.to.graphNEL(fNw)

}

