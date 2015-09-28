\name{keywordsTerm}
\alias{keywordsTerm}
\title{
Select keyword term
}
\description{
Selects a term as the most representative of the terms in the cluster based on keywords.
}
\usage{
keywordsTerm(termsDescriptions, nChar = 30)
}
\arguments{
  \item{termsDescriptions}{
List with the terms in each cluster. Output from \code{\link[FGNet:getTerms]{getTerms}}.
}
  \item{nChar}{
numeric. Maximum number of chars to show in the term. If the selected term is longer, it will be trimmed.
}
}
%\details{
%% Explain method?
%}
\value{
Character vector with the term selected for each cluster.
}
\seealso{
Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\examples{

# Previous Steps: FEA
results <- fea_gtLinker_getResults(3907019)

# Select keywords
keywordsTerm(getTerms(results), nChar=100)

}