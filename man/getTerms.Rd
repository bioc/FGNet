\name{getTerms}
\alias{getTerms}
\title{
Get terms in the metagroups/clusters.
}
\description{
Gets the terms in each metagroup/cluster (simplifyes the raw output from GeneTermLinker or DAVID).
}
\usage{
getTerms(feaResults, returnValue = "description")
}
\arguments{
  \item{feaResults}{
Output returned by any of the \link[FGNet:fea_david]{fea functions}.
}
  \item{returnValue}{
"description" Returns term description, e.g. "GO" returns GO term IDs.
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
results <- fea_gtLinker_getResults(jobID=1963186)
getTerms(results)
}
}
