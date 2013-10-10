\name{plotMetagroupsDistance}
\alias{plotMetagroupsDistance}
\title{
Plots distances between metagroups.
}
\description{
Plots the distances between metagroups taking into account the number of common genes.
}
\usage{
plotMetagroupsDistance(metagroupGenesMatrix)
}
\arguments{
  \item{metagroupGenesMatrix}{
  \emph{$metagroupGenesMatrix} returned by  \code{\link[FGNet:adjMatrix]{adjMatrix()}}.
}
}
\value{
 Plot and distance matrix.
}
\seealso{
	Full description of the package:
		\code{\link[FGNet:FGNet-package]{FGNet}}
}
\examples{
results <- getResults_gtLinker(jobID=1963186, jobName="gtLinker_ej")
adjMat <- adjMatrix(results$geneTermSets,
 attribute=results$metagroups[,"Silhouette Width", drop=FALSE])
plotMetagroupsDistance(adjMat$metagroupGenesMatrix)
}

