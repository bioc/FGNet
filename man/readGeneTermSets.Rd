\name{readGeneTermSets}
\alias{readGeneTermSets}
\title{
Read gene-term sets
}
\description{
Reads a file containing gene-term sets (formatted output from FEA) and transform it into clusters/metagroups and geneTermSets tables to imput to FGNet.
}
\usage{
readGeneTermSets(fileName, tool = NULL, simplifyGage = TRUE)
}
\arguments{
  \item{fileName}{
character. File name of the .txt file.
}
  \item{tool}{
character. Tool used for the FEA (row name from \code{data(FEA_tools); FEA_tools})
}
  \item{simplifyGage}{
logical. Only for gage: Wether to simplify the results and keep only the essential terms in the clustes.
}
}
\value{
List formated in the same way as the \link[FGNet:fea_david]{fea} functions. 
}
\seealso{
Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}

}