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
clustersDistance(incidenceMatices, mgCols = NULL)
}
\arguments{
  \item{incidenceMatices}{
  Object returned by  \code{\link[FGNet:fea2incidMat]{fea2incidMat()}}.
}
  \item{mgCols}{
  Colors for the metagroups.
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
results <- fea_gtLinker_getResults(jobID=1963186, jobName="gtLinker_example")
incidMatrices <- fea2incidMat(results)
clustersDistance(incidMatrices)
}
