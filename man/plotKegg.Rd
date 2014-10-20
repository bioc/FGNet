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

###################################
# Gene info
genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9", "CDC16", 
    "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1", "FIP1", 
    "GCN5", "GLC7", "HFI1", "KEM1", "LSM1", "LSM2", "LSM3", "LSM4", 
    "LSM5", "LSM6", "LSM7", "LSM8", "MPE1", "NGG1", "PAP1", "PAT1", 
    "PFS2", "PTA1", "PTI1", "REF2", "RNA14", "RPN1", "RPN10", "RPN11", 
    "RPN13", "RPN2", "RPN3", "RPN5", "RPN6", "RPN8", "RPT1", "RPT3", 
    "RPT6", "SGF11", "SGF29", "SGF73", "SPT20", "SPT3", "SPT7", "SPT8", 
    "TRA1", "YSH1", "YTH1")

# Gene expression (1=UP, -1=DW)
genesYeastExpr <- setNames(c(rep(1,29), rep(-1,30)), genesYeast)


###################################
# Plot pathway

# Specific pathway
plotKegg("sce04111", genesYeastExpr, geneIDtype="GENENAME")



# Pathways from analysis:
txtFile <- paste(file.path(system.file('examples', package='FGNet')), 
    "DAVID_Yeast_raw.txt", sep="/")
results <- format_david(txtFile, jobName="DavidAnalysis")

keggTerms <- unlist(getTerms(results, returnValue="KEGG"))

# plotKegg(keggTerms, ... geneIDtype="ENSEMBL")

}