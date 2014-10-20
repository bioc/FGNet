\name{plotGoAncestors}
\alias{plotGoAncestors}
\title{
Plot GO term ancestors
}
\description{
Plots the ancestors in the tree ontology for the given GO terms.
}
\usage{
plotGoAncestors(goIds, tColor = NULL, ontology = NULL,
    plotOutput = "static", nCharTerm = 50, nSize = NULL, labelCex = NULL, 
    asp = NULL, fileName = NULL, height = 1000)
}
\arguments{
  \item{goIds}{
character vector. GO IDs of the terms to plot.
}
  \item{tColor}{
character. Color for the term (i.e. based on expression).
}
  \item{ontology}{
character. If character determines which ontology to plot ("BP"", "MF" or "CC"").
}
  \item{plotOutput}{
"static", "dynamic" or "none". "static" will generate a standard R plot. "dynamic" will produce an interactive tkplot. "none" will not plot the network.
}
  \item{nCharTerm}{
numeric. Max term size (number of characters). Longer terms will be trimmed.
}
  \item{nSize}{
numeric. Determines the node size.
}
  \item{labelCex}{
numeric. Determines the node label size.
}
  \item{fileName}{
character. If provided, the plot is saved as png with this fileName.
}
  \item{asp}{
character. If fileName is provided, asp argument for plot.
}
  \item{height}{
numeric.  If fileName is provided, height argument for png().
}
}
\value{
An \link[base:invisible]{invisible} list with the nodes identified as leaves (\code{leaves}) and the graph (\code{iGraph}).
}

\seealso{
Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}

}
\examples{

plotGoAncestors(c("GO:0000152","GO:0043234", "GO:0044446", "GO:0043227"))

plotGoAncestors(c("GO:0051603", "GO:0019941", "GO:0051128","GO:0044265"), 
    plotOutput="dynamic")


# From analysis:
txtFile <- paste(file.path(system.file('examples', package='FGNet')),
    "DAVID_Yeast_raw.txt", sep="/")
results <- format_david(txtFile, jobName="DavidAnalysis")

plotGoAncestors(getTerms(results, returnValue="GO")$"Cluster 7", ontology="MF")
}