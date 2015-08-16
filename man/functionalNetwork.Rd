\name{functionalNetwork}
\alias{functionalNetwork}
\alias{intersectionNetwork}
\title{
Creates and plots the functional gene network.
}
\description{
Plots the functional networks. 

The default network links genes to genes, or terms to terms. The bipartite network links genes or terms to their clusters. 
}
\usage{
functionalNetwork(incidMatrices, plotType = c("default", "bipartite")[1], 
    plotOutput = "static", plotTitle = "Functional Network", 
    plotTitleSub = NULL, legendPrefix = NULL, legendText = NULL, 
    geneExpr = NULL, plotExpression = c("border", "fill"), 
    vExprColors=c(neg="#008000", zero="white", pos="#FF2020"),
    vSize = 12, vLabelCex = 2/3, vLayout = NULL, keepColors = TRUE, 
    bgTransparency = 0.4, eColor = "#323232", eWidth=NULL, weighted = FALSE, 
    keepAllNodes = FALSE, plotAllMg = FALSE)
}

\arguments{
  \item{incidMatrices}{
list or matrix. Raw output (list) from \code{\link[FGNet:fea2incidMat]{fea2incidMat}}: a list with slots:  "gtSetsMatrix", "filteredOut" and either "metagroupsMatrix" or "clustersMatrix". If only a matrix is provided, it will be asumed to be the clusters matrix, and all nodes will be connected to every other node in the metagroup.
}
  \item{plotType}{
"default" or "bipartite". 

Default network: Nodes are either genes or terms. Edges join nodes in common gene-term sets. Background and node color represent cluster/metagroup. White nodes are in several  clusters/metagroups.

Bipartite network: Nodes are genes or terms (circles) and their clusters (squares). 
By default it keeps only the genes or terms in more than one cluster or metagroup, which represents a simplified version of the functional network.
Node shape is only available in the "static" output.
}
  \item{plotOutput}{
"static", "dynamic" or "none". "static" will generate a standard R plot. "dynamic" will produce an interactive tkplot (metagroups background cannot be drawn). "none" will not plot the network.
}
  \item{plotTitle}{
character. Title to show on the plot.
}
  \item{plotTitleSub}{
character. Text to show at the bottom of the plot (sub-title).
}
  \item{legendPrefix}{
character. Label to show next to the cluster/metagroup id in the legend. In the bipartite network the legens replaces the cluster node label.
}
  \item{legendText}{
character. Description of each cluster (shown as the legend).  If FALSE, legend is not shown.
}
  \item{geneExpr}{
numeric. Named vector with the relative expression value of the gene (node). 0 is taken as reference, positive values will be plotted red, negative values green.
}
  \item{plotExpression}{
character. Determines the way to plot the expression: "border" adds a red or green border to the node, "fill" colors the whole with the expression color instead of the metagroup color.
}
  \item{vExprColors}{
character. Vector with the colors for expression: first color for negative values, second for zero, and third for positive.
}
  \item{vSize}{
numeric. Vertex size. If named, it allows to set a value for each gene. Name as "default" to set a default value, otherwise the default value is the mean.
}
  \item{vLabelCex}{
numeric. A numerical value giving the amount by which plotting text and symbols should be magnified relative to the default label size. 
}
  \item{vLayout}{
2 x n matrix or character. Where n is the number of nodes in the graph, each column gives the (x, y)-coordinates for the corresponding node. The bipartite network accepts "kk" (Kamada Kawai), "circle", or "sugiyama" (hierarquical).
}
  \item{keepColors}{
logical. If TRUE, it will keep the same colors for all the plots, independently of the filtered groups. Only available if incidMatrices is the raw result from \code{\link[FGNet:fea2incidMat]{fea2incidMat()}}.
}
  \item{bgTransparency}{
numeric. Value between 0 and 1 for the transparency of the metagroups background (only default network).
}
  \item{eColor}{
character. Color for the edges.
}
  \item{eWidth}{
numeric. Edge width. Not to plot edges, set eWidth=0 or eColor=NA.
}
  \item{weighted}{
logical. If TRUE, edges width will be based on the number of shared gene-term sets.
}
  \item{keepAllNodes}{
logical. Only used in bipartite network. If FALSE, nodes in only one cluster are not plotted. If TRUE, all nodes in the clusters are shown.
}
  \item{plotAllMg}{
logical. Only used in bipartite network. If FALSE, non-connected clusters are not plotted. If TRUE, all non-filtered clusters are shown.
}
}
\value{
Plots the functional networks.
 
An \link[base:invisible]{invisible} list with the igraph networks and incidence matrices, to collect it assign it to a variable.
}


\seealso{
Previous step in the workflow: \code{\link[FGNet:fea2incidMat]{fea2incidMat()}}
    
To see the terms included in each cluster or metagroup: \code{\link[FGNet:getTerms]{getTerms()}}

Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\examples{
###################################################
# Previous steps
# Set gene list:
genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9", 
    "CDC16", "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1", "FIP1", 
    "GCN5", "GLC7", "HFI1", "KEM1", "LSM1", "LSM2", "LSM3", "LSM4", "LSM5", 
    "LSM6", "LSM7", "LSM8", "MPE1", "NGG1", "PAP1", "PAT1", "PFS2", "PTA1", 
    "PTI1", "REF2", "RNA14", "RPN1", "RPN10", "RPN11", "RPN13", "RPN2", "RPN3", 
    "RPN5", "RPN6", "RPN8", "RPT1", "RPT3", "RPT6", "SGF11", "SGF29", "SGF73", 
    "SPT20", "SPT3", "SPT7", "SPT8", "TRA1", "YSH1", "YTH1")
# Optional gene expression
genesYeastExpr <- setNames(c(rep(1,29), rep(-1,30)), genesYeast) # 1=UP, -1=DW

# FEA:
# jobID <- query_gtLinker(genesYeast, organism = "Sc")
jobID <- 3907019
results <- fea_gtLinker_getResults(jobID)

###################################################
# Gene-based networks:
incidMat <- fea2incidMat(results, filterAttribute="Silhouette Width")

functionalNetwork(incidMat, geneExpr=genesYeastExpr)
functionalNetwork(incidMat, plotType="bipartite", plotOutput="dynamic", vSize=c(default=10, GLC7=20, PTA1=20))

getTerms(results)

# To modify the layout and plot as static network (with metagroup background)...
library(igraph)
# saveLayout <- tkplot.getcoords(1)   # tkp.id (ID of the tkplot window)
# functionalNetwork(incidMat, vLayout=saveLayout, plotType="bipartite")

# Only return the network, without plotting
fNw <- functionalNetwork(incidMat, plotOutput="none") 
class(fNw)
names(fNw)
betweenness(fNw$iGraph$commonClusters)

###################################################
# Term-based network
incidMat_terms <- fea2incidMat(results, key="Terms")
functionalNetwork(incidMat_terms, weighted=TRUE, plotOutput="dynamic")
functionalNetwork(incidMat_terms, plotType="bipartite", plotOutput="dynamic", 
    plotAllMg=TRUE)
functionalNetwork(incidMat_terms, plotType="bipartite", plotOutput="dynamic",
    keepAllNodes=TRUE)

# Including generic terms filterd by GtLinker from final metagroups:
incidMat_terms2 <- fea2incidMat(results, key="Terms", removeFiltered=FALSE)
functionalNetwork(incidMat_terms2, weighted=TRUE)
}
