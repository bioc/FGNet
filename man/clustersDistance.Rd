\name{clustersDistance}
\alias{clustersDistance}
\alias{plotMetagroupsDistance}
\title{
Plots distances between metagroups.
}
\description{
Plots the distances between metagroups taking into account the number of common genes.
}
\usage{
clustersDistance(incidenceMatices, mgCols = NULL, clustMethod="average")
}
\arguments{
  \item{incidenceMatices}{
  Object returned by  \code{\link[FGNet:fea2incidMat]{fea2incidMat()}}.
}
  \item{mgCols}{
  Colors for the metagroups.
}
  \item{clustMethod}{
  Clustering method. Character string (i.e. "single", "complete", "average") for function \code{\link[stats:hclust]{hclust()}} (argument 'method').
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
\dontrun{
results <- fea_gtLinker_getResults(jobID=1963186, jobName="gtLinker_example")
incidMatrices <- fea2incidMat(results)
clustersDistance(incidMatrices)
}
}
