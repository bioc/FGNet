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

# Previous Steps
genesYeast <- c("YBL084C", "YDL008W", "YDR118W", "YDR301W", "YDR448W", "YFR036W", 
    "YGL240W", "YHR166C", "YKL022C", "YLR102C", "YLR115W", "YLR127C", "YNL172W", 
    "YOL149W", "YOR249C")
results <- fea_david(genesYeast) # email=example@email.com

# Select keywords
keywordsTerm(getTerms(results), nChar=100)

}