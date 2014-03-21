\name{functionalNetwork}
\alias{functionalNetwork}
\title{
Create and plot the functional gene network.
}
\description{
Plots the gene network. Edges join nodes with  common gene-term sets. Background and node color represent metagroup/clusters. White nodes are in several metagroups/clusters.
}
\usage{

functionalNetwork(incidMatrices, plotType="static", returnGraph=FALSE,
plotTitle="Functional Network", vSize=12, vLabelCex=3/4, 
vLayout=NULL, bgTransparency=0.4, legendMg=NULL, keepColors=TRUE,
eColor="#323232", weighted=FALSE)
}

\arguments{
 \item{incidMatrices}{
list or matrix. Raw output (list) from \code{\link[FGNet:toMatrix]{toMatrix}}.  List with slots:  "gtSetsMatrix", "filteredOut" and either "metagroupsMatrix" or "clustersMatrix". If only a matrix is provided, it will be asumed to be the metagroups matrix, and all nodes will be connected to every other node in the metagroup.
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
  \item{eColor}{
character. Color for the edges.
}
  \item{weighted}{
logical. If TRUE, edges width will be based on the number of shared gene-term sets.
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
incidMat <- toMatrix(results, threshold=0,
		attribute=results$metagroups[,"Silhouette Width", drop=FALSE])

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


# Term-based network
incidMatTerms <- toMatrix(results, key="Terms")
functionalNetwork(incidMatTerms, plotType="dynamic", weighted=TRUE, eColor="grey")

# Including generic terms filterd by GtLinker from final metagroups:
incidMatTerms2 <- toMatrix(results, key="Terms", removeFiltered=FALSE)
functionalNetwork(incidMatTerms2, weighted=TRUE)

### DAVID:
genesMetabolism <- c("YGR175C", "YHR007C", "YMR202W", "YJL167W", 
 "YNL280C", "YGR060W", "YGL001C", "YLR100W", "YLR056W", "YGL012W",
 "YMR015C", "YML008C", "YHR072W", "YHR190W", "YKL004W", "YBR036C",
 "YDR294C", "YDR072C", "YKL008C", "YHL003C", "YMR296C", "YDR062W",
 "YJL134W", "YOR171C", "YLR260W", "YMR298W", "YMR272C", "YPL057C",
 "YDR297W", "YBR265W", "YPL087W", "YBR183W", "YKR053C")
txtFile <- query_david(genesMetabolism)
results <- getResults_david(txtFile)
														
# Term-based network
incidMatTerms <- toMatrix(results, key="Terms")
functionalNetwork(incidMatTerms$clustersMatrix, plotType="dynamic")
									


}
