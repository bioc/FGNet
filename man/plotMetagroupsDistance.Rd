\name{plotMetagroupsDistance}
\alias{plotMetagroupsDistance}
\title{
Plots distances between metagroups.
}
\description{
Plots the distances between metagroups taking into account the number of common genes.
}
\usage{
plotMetagroupsDistance(incidenceMatices)
}
\arguments{
  \item{incidenceMatices}{
  Object returned by  \code{\link[FGNet:toMatrix]{toMatrix()}}.
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
incidMat <- toMatrix(results,
 attribute=results$metagroups[,"Silhouette Width", drop=FALSE])
plotMetagroupsDistance(incidMat)
}
