\name{getResults_gtLinker}
\alias{getResults_gtLinker}

\title{
Get results from GeneTerm Linker.
}
\description{
Retrieves the metagroups and gene-term sets for the given jobID from GeneTerm Linker. The jobID can be obtained from \link[FGNet:query_gtLinker]{query_gtLinker} or from an analysis performed directly at GeneTerm Linker web (http://gtlinker.cnb.csic.es) . 
}
\usage{
getResults_gtLinker(jobID, path = getwd(), jobName = "", 
alreadyDownloaded = FALSE,  keepTrying=FALSE,
serverWeb = "http://gtlinker.cnb.csic.es",
serverWS =  "http://gtlinker.cnb.csic.es:8182")
}
\arguments{
  \item{jobID}{ 
  numeric. ID of the job/analysis in GeneTerm Linker.
}
  \item{path}{
  character. Directory in which to save the files.
}
  \item{jobName}{
  character. Folder name and prefix for the files.
}
  \item{alreadyDownloaded}{
  logical. If the files have already been downloaded, these will be read. Make sure to use the appropiate jobID and job name.
}
  \item{keepTrying}{
  logical. If true, if the job has not finished, it will keep trying to get the results every few seconds.
}
  \item{serverWeb}{
character. GeneTerm Linker web server. It should match the web service or web address in which the analysis was performed. \cr
Available mirrors: "http://gtlinker.cnb.csic.es" and "http://cicblade.dep.usal.es:8000"
}
  \item{serverWS}{
character. GeneTerm Linker webservice server. \cr
Available mirrors: "http://gtlinker.cnb.csic.es:8182" and "http://cicblade.dep.usal.es:8182". If you change the webserice server, make sure to use the matching 'serverWeb'.
}

}
\value{
List:\cr
\code{globalMetagroups } data.frame containing the metagroups and their information:
	\itemize{ 
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
}
\seealso{
	Next step in the workflow: 
		\code{\link[FGNet:adjMatrix]{adjMatrix()}}\cr
		
	Previous step in the workflow: 
		\code{\link[FGNet:query_gtLinker]{query_gtLinker()}}\cr
	
	Equivalent function for DAVID:
	\code{\link[FGNet:getResults_david]{getResults_david}}\cr
	
	Full description of the package:
	\code{\link[FGNet:FGNet-package]{FGNet}}
}

\examples{
jobID <- 3907019
results <- getResults_gtLinker(jobID)
}
