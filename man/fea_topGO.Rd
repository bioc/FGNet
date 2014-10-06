\name{fea_topGO}
\alias{fea_topGO}
\title{
FEA - topGO
}
\description{
Performs the functional enrichment analysis through topGO [1].
}
\usage{
fea_topGO(geneList, geneIdType = "ENSEMBL", geneLabels=NULL, organism = "Hs", 
    annotations = c("GO_BP", "GO_MF", "GO_CC"), 
    genesUniverse = NULL, refPackage = NULL, 
    geneID2GO = NULL, nodeSize = 5, pValThr = 0.01, testStat = NULL, 
    jobName = NULL)
}
\arguments{
  \item{geneList}{
character vector. List of genes to analyze.
}
  \item{geneIdType}{
character. Type of gene identifier should be available for the organism package.
}
\item{geneLabels}{
    named character vector. Gene name or label to use in the report/plots instead of the original gene ID. The vector names should be the gene ID and the content of the vector the gene label. The resulting geneTermSets table will contain the original gene ID column (geneIDs) and the label column (Genes).
}
  \item{organism}{
two letter code for the organism. See: \code{data(organisms);organisms}
}
  \item{annotations}{
character vector. Annotation spaces for the functional analysis. Accepted values: "GO_BP", "GO_MF", "GO_CC".
}
  \item{genesUniverse}{
character vector. List of genes used for background (i.e. all genes available in the chip).
}
  \item{refPackage}{
character. Name of the package to use for calculating the genes universe. A Chip package is recommended. If NULL the genes universe is set as all the genes available in the organism package.
}
  \item{geneID2GO}{
GO gene sets. If NULL it is calculated automatically. %To speed up the execution of multiple queries you can use: \code{geneID2GO <- buildDatabases(orgPackage, geneIdType)$geneID2GO}
}
  \item{nodeSize}{
numeric. Minimum size of GO terms. TopGo authors recommend 5-10 for more stable results, 1 for no prune.
}
  \item{pValThr}{
numeric. P-value threshold.
}
  \item{testStat}{
classicCount from toGO. If NULL: GOFisherTest is used.
}
  \item{jobName}{
character. Folder name and prefix for the files.
}
}

\value{
\link[base:invisible]{Invisible} list with the folowing fields:\cr
\code{queryArgs } list with the arguments for the query.

\code{clusters } Empty list. only for compatibility.
 
\code{geneTermSets } data.frame containing the gene-term sets.
\itemize{ 
    \item Ont: Ontology to wich the term belongs (BP, MF or CC)
    \item Terms: Term in the gene-term set.
    \item Genes: Genes in the gene-term set.
    \item GenesIDs: In case GeneLabels was provided, original gene ID.
    \item Other stats provided by \link[topGO:getSigGroups]{topGO}: Annotated, Significant, Expected, classic.
}
 
\code{fileName} .txt file with the formatted FEA results.
}
\seealso{
Other FEA tools:
\itemize{
    \item{\code{\link[FGNet:fea_david]{fea_david()}} (Requires internet connection)}
    \item{\code{\link[FGNet:fea_gtLinker]{fea_gtLinker()}} & \code{fea_gtLinker_getResults()}
          (Requires internet connection)}
    \item{\code{\link[FGNet:fea_gage]{fea_gage()}}
    
   To import results from a previous/external FEA analysis: \code{\link[FGNet:format_david]{format_david()}}, \code{format_results()} and \code{\link[FGNet:readGeneTermSets]{readGeneTermSets()}}.}
}

Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}

\references{
[1] Adrian Alexa and Jorg Rahnenfuhrer (2010) topGO: Enrichment analysis for Gene Ontology. R package version 2.16.0. URL: \url{http://www.bioconductor.org/packages/release/bioc/html/topGO.html}
}

\examples{
\dontrun{
    

# Load/format gene list:
genesYeast <- c("YBL084C", "YDL008W", "YDR118W", "YDR301W", "YDR448W", "YFR036W",
    "YGL240W", "YHR166C", "YKL022C", "YLR102C", "YLR115W", "YLR127C", "YNL172W",
    "YOL149W", "YOR249C")
geneExpr <- setNames(c(rep(1,10),rep(-1,5)), genesYeast) 

library(org.Sc.sgd.db)
geneLabels <- unlist(as.list(org.Sc.sgdGENENAME)[genesYeast])
names(geneExpr) <- geneLabels[names(geneExpr)] 

# FEA:
results <- fea_topGO(genesYeast, geneIdType="ENSEMBL", 
    geneLabels=geneLabels, organism="Sc") 

# To continue the workflow... (see help for further details)
FGNet_report(results, geneExpr=geneExpr)

incidMat <- fea2incidMat(results)
functionalNetwork(incidMat)
    
    
}
}