\name{fea_david}
\alias{fea_david}
\alias{query_david}
\title{
FEA - DAVID
}
\description{
Performs the functional enrichment analysis and clustering through DAVID [1] (requires internet connection).
}
\usage{
fea_david(geneList, geneIdType = "ENSEMBL_GENE_ID", geneLabels=NULL, 
  annotations = c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", 
 "INTERPRO"), 
  email = NULL, 
  argsWS = c(overlap = 4L, initialSeed = 4L, finalSeed = 4L, linkage = 0.5, 
  kappa = 35L), jobName = NULL, downloadFile=TRUE)
}
\arguments{
  \item{geneList}{
character vector. List of genes to analyze.
}
  \item{geneIdType}{
character vector. Type of gene identifier.\cr
Web API: ENSEMBL_GENE_ID, ENTREZ_GENE_ID, GENE_SYMBOL, UNIPROT_ID... For more, check DAVID's API documentation.\cr
Web Service: run \code{getIdTypes(DAVIDWebService$new(email=...))}
}
\item{geneLabels}{
    named character vector. Gene name or label to use in the report/plots instead of the original gene ID. The vector names should be the gene ID and the content of the vector the gene label. The resulting geneTermSets table will contain the original gene ID column (geneIDs) and the label column (Genes).
}
  \item{annotations}{
character vector. Annotation spaces for the functional analysis. \cr
Web API: check DAVID's API documentation.\cr
Web Service: run \cr
\code{getAllAnnotationCategoryNames(DAVIDWebService$new(email=...))}.
}
  \item{email}{
character. If provided, the query will be performed though DAVID's Web Service (recommended). Requires registration (see details).
}
  \item{argsWS}{
named integer vector. Additional arguments for the clustering. Only available using the web service. 
}
  \item{jobName}{
character. Folder name and prefix for the files.
}
  \item{downloadFile}{
logical. If TRUE, the result files are saved in the current directory (required to generate report).
}
}
\details{
To perform the queries, please register at \url{http://david.abcc.ncifcrf.gov/webservice/register.htm}.\cr

\bold{NOTE:} Since August 2015, DAVID requires https. This causes errors in some systems. 
A (hopefully) temporary solution requries to install some certificates locally. See RDAVIDWebService help: \url{https://support.bioconductor.org/p/70090/#72226}\cr


As an alternative, the web API allows to perform a small query without registering. \cr
Note this option is not available in some systems, and the maximum number of genes is limited to 400. (It can be less depending on the ID types and the length of the resulting URL).\cr
More details and full list of gene ID types and annotations are available at: 
\url{http://david.abcc.ncifcrf.gov/content.jsp?file=DAVID_API.html}.\cr

If the functional annotation and clustering has been performed directly at DAVID's \bold{Website} (\url{http://david.abcc.ncifcrf.gov/summary.jsp}) fea_david() is not required. Instead, provide the file (or the URL of the file) containing the results of the analysis to \code{\link[FGNet:format_david]{format_david()}}.\cr
}



\section{Warning}{
The web service and the API have different default arguments. To obtain the same results with both methods use:

\code{API_defaults <- c(overlap=3L, initialSeed=3L, finalSeed=3L, linkage=0.5, kappa=50L)}

\code{fea_david(genesYeast, email="example@email.com", argsWS=API_defaults)}
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
    \item Category: Type of annotation of the term (i.e. GO, ...)
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
Other FEA tools:
\itemize{
    \item{\code{\link[FGNet:fea_gtLinker]{fea_gtLinker()}} & \code{fea_gtLinker_getResults()}
          (Requires internet connection)}
    \item{\code{\link[FGNet:fea_gage]{fea_gage()}}}
    \item{\code{\link[FGNet:fea_topGO]{fea_topGO()}} 
    
   To import results from a previous/external FEA analysis: \code{\link[FGNet:format_david]{format_david()}}, \code{format_results()} and \code{\link[FGNet:readGeneTermSets]{readGeneTermSets()}}.}
}

Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\examples{


# Load/format gene list:
geneList <- c("YBL084C", "YDL008W", "YDR118W", "YDR301W", "YDR448W", "YFR036W",
    "YGL240W", "YHR166C", "YKL022C", "YLR102C", "YLR115W", "YLR127C", "YNL172W", 
    "YOL149W", "YOR249C")

library(org.Sc.sgd.db)
geneLabels <- unlist(as.list(org.Sc.sgdGENENAME)[geneList])

geneExpr <- setNames(c(rep(1,10),rep(-1,5)), geneLabels) 

\dontrun{
# DAVID
results_David <- fea_david(geneList, geneLabels=geneLabels, email="example@email.com")

# Available IDs and annotations:
getIdTypes(DAVIDWebService$new(email="example@email.com"))
getAllAnnotationCategoryNames(DAVIDWebService$new(email="example@email.com"))

results <- fea_david(geneList, geneIdType="ENSEMBL_GENE_ID",
    annotations="GOTERM_BP_ALL", email="example@email.com", jobName="yeastDavid")


# To continue the workflow... (see help for further details)
FGNet_report(results, geneExpr=geneExpr)

incidMat <- fea2incidMat(results)
functionalNetwork(incidMat, geneExpr=geneExpr)


}
}
