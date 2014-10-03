\name{analyzeNetwork}
\alias{analyzeNetwork}
\title{
Analyze Functional Network
}
\description{
Analyzes the degree and betweeness of the genes in the functional network.
}
\usage{
analyzeNetwork(incidMatrices, fNw = NULL, plotOutput = TRUE, colors = NULL)
}
\arguments{
  \item{incidMatrices}{
list or matrix. Output from \code{\link[FGNet:fea2incidMat]{fea2incidMat}} or the equivalent incidence matrices.
}
  \item{fNw}{
list. Return from \code{\link[FGNet:functionalNetwork]{functionalNetwork}} to avoid recalculating.
}
  \item{plotOutput}{
logical. Wether to plot the degree and betweeness boxplots.
}
  \item{colors}{
vector. Colors for the metagroups
}
}
\value{
List:

\itemize{ 
    \item degree, betweeness: Degree and Betweeness of the nodes in the global network (commonClusters) and within each cluster/metagroup (subsets of commonGtSets network).
    \item transitivity: Transitivity of the networks.
    \item hubsTable: Betweeness of each node in each cluster.
    \item hubsList: Nodes selected as hubs in the global network and within each cluster/metagroup (nodes with betweeness over 75\% in the given network/subnetwork).
    \item intraHubsCount: Number of times each node was selected as intra-cluster hub.
}
}

\seealso{
Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}

}
\examples{
# Previous Steps
jobID <- 3907019
results <- fea_gtLinker_getResults(jobID)
incidMat <- fea2incidMat(results, filterAttribute="Silhouette Width")


# Plot node degree and betweensess
analyzeNetwork(incidMat)

# Get stats without plotting
nwStats <- analyzeNetwork(incidMat, plotOutput=FALSE)
names(nwStats)
nwStats$hubsTable

}
