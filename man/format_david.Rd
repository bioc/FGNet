\name{format_david}
\alias{format_david}
\alias{getResults_david}
\title{
Format DAVID output
}
\description{
Format DAVID 'functional annotation and clustering' output to use with FGNet.
}
\usage{
format_david(fileName, 
jobName = NULL, 
geneLabels = NULL, 
moveFile = FALSE, 
downloadFile=TRUE)
}
\arguments{
  \item{fileName}{
character. URL or local file with the results of a DAVID analysis, for example, performed at DAVID's Website (\url{http://david.abcc.ncifcrf.gov/summary.jsp}).
In case of local file, it should be the absolute path to the .txt file (whole location from root to the file including file name: "C:\\\\Documents\\\\23424203.txt", "/home/user/2342342.txt").
}
  \item{jobName}{
character. Folder name and prefix for the formatted files.
}
\item{geneLabels}{
named character vector. Gene name or label to use in the report/plots instead of the original gene ID. The vector names should be the gene ID and the content of the vector the gene label. The resulting geneTermSets table will contain the original gene ID column (geneIDs) and the label column (Genes).
}
\item{moveFile}{  
logical. If TRUE the original file is moved to the new location. If FALSE, the file is copied.
}
  \item{downloadFile}{
logical. If TRUE, the result files are saved in the current directory (required to generate report).
}
}
\value{
\link[base:invisible]{Invisible} list with the folowing fields:\cr
\code{queryArgs } list with the arguments for the query.

\code{clusters } data.frame containing the clusters and their information:
\itemize{ 
    \item Cluster: Cluster ID.
    \item nGenes: Number of genes in the cluster.
    \item ClusterEnrichmentScore: Score for the cluster.
    \item Genes: Genes in the cluster.
    \item Terms: Terms in the cluster.
    \item keyWordsTerm: Term is the most representative of the terms in the cluster based on keywords.
}
 
\code{geneTermSets } data.frame containing the gene-term sets that support each cluster.
\itemize{ 
    \item Cluster: Number (id) of the cluster the gene-term set belongs to.
    \item ClusterEnrichmentScore: Score for the cluster. Same value for all terms in each cluster.
    \item Category: Type of annotation of the term (i.e. GO, Kegg...)
    \item Terms: Term in the gene-term set.
    \item Genes: Genes in the gene-term set.
    \item GenesIDs: In case GeneLabels was provided, original gene ID.
    \item Other stats: Count, PValue, List.Total, Pop.Hits, Pop.Total, Fold.Enrichment, Bonferroni, Benjamini, FDR.
}
 
\code{fileName} .txt file with the formatted FEA results.
}
\references{
[1] Huang DW, Sherman BT, Lempicki RA (2009) Bioinformatics enrichment tools: paths toward the comprehensive functional analysis of large gene lists. Nucleic Acids Res. 37(1):1-13.
}
\seealso{
Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\examples{

# Select file:
txtFile <- "http://david.abcc.ncifcrf.gov/data/download/901234901248.txt"
txtFile <- paste(file.path(system.file('examples', package='FGNet')), 
    "DAVID_Yeast_raw.txt", sep="/")

# Read:
results <- format_david(txtFile, jobName="DavidAnalysis")


# To continue the workflow... (see help for further details)
getTerms(results)
incidMat <- fea2incidMat(results)
functionalNetwork(incidMat)

?FGNet_report

}
