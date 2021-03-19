\name{fea_gtLinker}
\alias{fea_gtLinker}
\alias{fea_gtLinker_getResults}
\alias{query_gtLinker}
\alias{getResults_gtLinker}
\title{
FEA - Gene-Term Linker
}
\description{
Performs the functional enrichment analysis and clustering through Gene-Term Linker [1] (requires internet connection).

Since Gene-Term Linker takes a while to analyze the gene list, the process has been splitted in two steps:
\enumerate{
\item{fea_gtLinker(): Submits the query}
\item{fea_gtLinker_getResults(): Retrieves the results of the analysis. It might take a few minutes for the results to become available.}
}
}
\usage{
fea_gtLinker(geneList, organism = "Hs", 
    annotations = c("GO_Biological_Process", "GO_Molecular_Function", 
    "GO_Cellular_Component", "InterPro_Motifs"),
    minSupport = 4, serverWS = "http://gtlinker.cnb.csic.es:8182")

fea_gtLinker_getResults(jobID = NULL, organism = NULL, jobName = NULL, 
    alreadyDownloaded = FALSE, keepTrying = FALSE, 
    serverWeb = "http://gtlinker.cnb.csic.es", 
    serverWS = "http://gtlinker.cnb.csic.es:8182")
}
\arguments{

\bold{fea_gtLinker():}

  \item{geneList}{
character vector. List of genes to analyze.
}
  \item{annotations}{
character vector. Annotation spaces for the functional analysis. \cr
Available values: "GO_Biological_Process", "GO_Molecular_Function", \cr
"GO_Cellular_Component", "InterPro_Motifs".
}
  \item{minSupport}{
numeric. Minimum number of genes per group.
}
\bold{common arguments}:
  \item{serverWS}{
character. GeneTerm Linker webservice server. \cr
Available mirrors: "http://gtlinker.cnb.csic.es:8182" %and "http://cicblade.dep.usal.es:8182".
If you change the webserice server, make sure to use the matching 'serverWeb' in the following step.
}
  \item{organism}{
character. "Hs" (Homo sapiens) or "Sc" (Saccharomyces cerevisiae).
}

\bold{fea_gtLinker_getResults()}:

  \item{jobID}{
numeric. ID of the job/analysis in GeneTerm Linker.
}
  \item{jobName}{
character. Folder name and prefix for the files.
}
  \item{alreadyDownloaded}{
logical. If the files have already been downloaded, these will be read instead of downloaded again.
}
  \item{keepTrying}{
logical. If true, if the job has not finished, it will keep trying to get the results every few seconds.
}
  \item{serverWeb}{
character. GeneTerm Linker web server. It should match the web service or web address in which the analysis was performed. \cr
Available mirrors: "http://gtlinker.cnb.csic.es" % and "http://cicblade.dep.usal.es:8000"
}

}
\value{

\bold{fea_gtLinker()} returns the jobID of the analysis

\bold{fea_gtLinker_getResults()} returns an \link[base:invisible]{invisible} list with the folowing fields:\cr

\code{queryArgs } list with the arguments for the query.

\code{metagroups } data.frame containing the metagroups and their information:
\itemize{ 
    \item Metagroup: Metagroup ID.
    \item Size: Number of gene-term sets supporting the metagroup.
    \item Diameter: Maximum Cosine distance within the GeneTerm-sets of each metagroup (ranges from 0 to 1).
    \item Similarity: 1 - average Cosine distance within the GeneTerm-sets of each metagroup (ranges from 0 to 1). Distance and similarity calculations are done based on the genes present in the metagroups.
    \item Silhouette Width: Measures the compactness and proximity of multiple groups (ranges from 1 to -1). 
Metagroups with negative Silhouette Width usually include diverse annotations and genes with low functional coherence. 
    \item Genes: Genes in the metagroup.
    \item nGenes: Number of genes in the metagroup.
    \item nref_list: Number of annotated genes in the reference list.
    \item pValue: Adjusted p-value.
    \item Terms: Non-generic terms in the metagroup.
 }
\code{geneTermSets }data.frame containing the gene-term sets that support each metagroup.
	\itemize{ 
	\item Metagroup: Id of the metagroup the gene-term set belongs to.
	\item Genes: Genes in the gene-term set.
	\item nGenes: Number of annotated genes in the input list. In brackets: Total number of genes in the input list.
	\item nref_list: Number of annotated genes in the reference list.  In brackets: Total number of genes in the reference list.
	\item pValue: Adjusted p-value.
	\item Terms: Terms in the gene-term set.
 }
 
\code{fileName} .txt file with the formatted FEA results.
}
\references{
[1] Fontanillo C, Nogales-Cadenas R, Pascual-Montano A, De Las Rivas J (2011) Functional Analysis beyond Enrichment: Non-Redundant Reciprocal Linkage of Genes and Biological Terms. PLoS ONE 6(9): e24289. URL: \url{http://gtlinker.cnb.csic.es}
}

\seealso{
Other FEA tools:
\itemize{
    \item{\code{\link[FGNet:fea_gage]{fea_gage()}}}
    \item{\code{\link[FGNet:fea_topGO]{fea_topGO()}} 
    
   To import results from a previous/external FEA analysis: \code{\link[FGNet:format_david]{format_david()}}, \code{format_results()} and \code{\link[FGNet:readGeneTermSets]{readGeneTermSets()}}.}
}

Overview of the package: \code{\link[FGNet:FGNet-package]{FGNet}}

Package tutorial: \code{vignette("FGNet-vignette")}
}
\examples{
\dontrun{
### Execute FEA:
genesYeast <- c("ADA2", "APC1", "APC11", "APC2", "APC4", "APC5", "APC9", 
     "CDC16", "CDC23", "CDC26", "CDC27", "CFT1", "CFT2", "DCP1", "DOC1", "FIP1", 
     "GCN5", "GLC7", "HFI1", "KEM1", "LSM1", "LSM2", "LSM3", "LSM4", "LSM5", 
     "LSM6", "LSM7", "LSM8", "MPE1", "NGG1", "PAP1", "PAT1", "PFS2", "PTA1", 
     "PTI1", "REF2", "RNA14", "RPN1", "RPN10", "RPN11", "RPN13", "RPN2", "RPN3", 
     "RPN5", "RPN6", "RPN8", "RPT1", "RPT3", "RPT6", "SGF11", "SGF29", "SGF73", 
     "SPT20", "SPT3", "SPT7", "SPT8", "TRA1", "YSH1", "YTH1")
# Optional expression (1=UP, -1=DW):
genesYeastExpr <- setNames(c(rep(1,29), rep(-1,30)), genesYeast)

# Submit query
jobID <- fea_gtLinker(geneList=genesYeast,organism="Sc")
jobID


### Get results from FEA:
jobID <- 3907019 # job ID of the query
results <- fea_gtLinker_getResults(jobID=jobID)

# To continue the workflow... (see help for further details))
incidMat <- fea2incidMat(results)
functionalNetwork(incidMat)

# Or full report
FGNet_report(results, geneExpr=genesYeastExpr)
}

}
