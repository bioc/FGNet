\name{plotKegg}
\alias{plotKegg}
\title{
Plot KEGG pathway
}
\description{
Plots KEGG pathway with the given genes.
}
\usage{
plotKegg(keggIDs, geneExpr, geneIDtype = "ENSEMBL", 
    colType = c("continuous", "discrete"))
}
\arguments{
  \item{keggIDs}{
character vector. KEGG IDS with prefix.
}
  \item{geneExpr}{
named numeric vector. Names should contain the gene ID, the value the gene expression. NA is accepted as expression value.
}
  \item{geneIDtype}{
character vector. Type of gene identifier. 
}
  \item{colType}{
"continuous" or "discrete" determines the color of the genes.
}
}
\details{
Color code:
\itemize{
\item Green: Negative expression.
\item Red: Positive expression.
\item Light yellow: Expression zero or very close to zero.
\item Blue: Expression not available.
\item White: Gene not in list.
}
}
\value{
The pathway plot is saved in the current working directory.
}
\seealso{
Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}

}
\examples{

geneList<- c("YPL087W", "YDR297W", "YBR183W", "YMR296C", "YKL004W", "YPL057C", 
    "YOR171C", "YMR272C", "YHL003C", "YDR072C", "YMR298W", "YBR036C", "YLR260W", 
    "YKL008C", "YDR294C", "YBR265W", "YJL134W", "YKR053C", "YDR062W", "YML008C", 
    "YGL001C", "YHR190W", "YGL012W", "YMR015C", "YJL167W", "YHR072W", "YLR100W", 
    "YMR202W", "YLR056W", "YGR060W", "YGR175C", "YHR007C", "YNL280C")
geneExpr <- setNames(sample(c(-1,0,1,NA), length(geneList), replace=TRUE), geneList) 

plotKegg("sce00100", geneExpr, geneIDtype="ENSEMBL")


# From analysis:
txtFile <- paste(file.path(system.file('examples', package='FGNet')), 
    "David_Metabolism.txt", sep="/")
results <- format_david(txtFile, jobName="DavidAnalysis")

keggTerms <- unlist(getTerms(results, returnValue="KEGG"))

plotKegg(keggTerms, geneExpr, geneIDtype="ENSEMBL")

}