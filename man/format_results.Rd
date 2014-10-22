\name{format_results}
\alias{format_results}
\title{
Format FEA results from external tools. 
}
\description{
Format the functional analysis results from external tools to use with FGNet. 
}
\usage{
format_results(fileName, newFileName = NULL, clusterCol = NULL, 
    geneCol = NULL, geneSep = NULL, termDescCol = NULL, termIDCol = NULL,
    termCatCol = NULL, termCat = NULL, termSep = NULL, 
    tool = "Imported text file", simplifyGage = TRUE, ...)
}
\arguments{
  \item{fileName}{
character. File name with the FEA results.
}
  \item{newFileName}{
character. Name for the formatted files.
}
  \item{clusterCol}{
character. Name of the column to use for clustering.
}
  \item{geneCol}{
character. Name of the column with the genes.
}
  \item{geneSep}{
character. Character separating diferent genes in the same field (i.e. ",", ";", ...)
}
  \item{termDescCol}{
character. Name of the column with the terms description.
}
  \item{termIDCol}{
character. Name of the column with the terms ID.
}
  \item{termCatCol}{
character. Name of the column with the terms type/category.
}
  \item{termCat}{
character. Name of the annotation type if it is common to all gene-term sets. Provide either termCatCol or termCat, not both.
}
  \item{termSep}{
character. Character separating diferent terms in the same field (i.e. ",", ";", ...)
}
  \item{tool}{
character. Tool used for the FEA (row name from \code{data(FEA_tools); FEA_tools})
}
  \item{simplifyGage}{
logical. For internal use, only for GAGE. Determines wether to keep non essential terms in the final clusters.
}
  \item{\dots}{
Further argumets to pass to "read.table"
}
}
\value{
Saves the formatted file and returns an  \link[base:invisible]{invisible} list with the appropiate format to use with \link[FGNet:FGNet_report]{FGNet_report()} and \link[FGNet:fea2incidMat]{fea2incidMat()} (fields "clusters", "geneTermSets" and "fileName").
}
\seealso{
Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\examples{

\dontrun{

results <- format_results("/home/user/feaResults.txt", clusterCol="Cluster", 
    geneCol="Genes", termDescCol="Terms", sep="\t")

}
}