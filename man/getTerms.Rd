\name{getTerms}
\alias{getTerms}
\title{
Get terms in the metagroups/clusters.
}
\description{
Gets the terms in each metagroup/cluster (simplifyes the raw output from GeneTermLinker or DAVID).
}
\usage{
getTerms(results)
}
\arguments{
  \item{results}{
Output returned by \code{\link[FGNet:getResults_gtLinker]{getResults_gtLinker()}} or \code{\link[FGNet:getResults_david]{getResults_david()}}.
}
}
\value{
\item{List of matrices}{Each matrix contais the terms in each metagroup. This matrix contains only the term description. To get the term ID, annotation type, number of genes, or any other information, see the raw results returned by getResults.}
}
\seealso{		
Full description of the package:
\code{\link[FGNet:FGNet-package]{FGNet}}
}
\examples{
\dontrun{
# GeneTermLinker
jobID <- 9396814
results <- getResults_gtLinker(jobID)
getTerms(results)
}
}
